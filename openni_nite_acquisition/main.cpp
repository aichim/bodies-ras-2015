#include <array>
#include <iostream>
#include <map>
#include <vector>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <OpenNI.h>
#include <NiTE.h>
#include <fstream>


using namespace std;
using namespace openni;
using namespace nite;


void
getJointPos (const nite::Skeleton &skeleton,
			 std::vector<cv::Vec3f> &joint_pos,
			 std::vector<float> &confidences)
{
	std::vector<nite::JointType> types = {JOINT_TORSO, JOINT_NECK, JOINT_LEFT_SHOULDER, JOINT_RIGHT_SHOULDER, JOINT_LEFT_ELBOW, JOINT_RIGHT_ELBOW,	
										  JOINT_LEFT_HAND, JOINT_RIGHT_HAND, JOINT_HEAD, JOINT_LEFT_HIP, JOINT_RIGHT_HIP, JOINT_LEFT_KNEE, JOINT_RIGHT_KNEE,
										  JOINT_LEFT_FOOT, JOINT_RIGHT_FOOT};
    for (size_t j_i = 0; j_i < types.size (); ++j_i)
    {
		confidences.push_back (skeleton.getJoint (types[j_i]).getPositionConfidence ());
		cv::Vec3f pos;
		pos[0] = skeleton.getJoint (types[j_i]).getPosition ().x;
		pos[1] = skeleton.getJoint (types[j_i]).getPosition ().y;
		pos[2] = skeleton.getJoint (types[j_i]).getPosition ().z;
		joint_pos.push_back (pos);
	}
}


void
visualizeSkeleton (const std::vector<cv::Vec3f> &joint_pos, 
				   cv::Mat &image)
{
	for (size_t i = 0; i < joint_pos.size (); ++i)
	{
		double u = 525. * joint_pos[i][0] / joint_pos[i][2] + 319.5;
		double v = 480. - (525. * joint_pos[i][1] / joint_pos[i][2] + 239.5);

		cv::circle (image, cv::Point (u, v), 5, cv::Scalar (255, 255, 255));
	}
}


void
saveSkeletonFile (std::string str, 
				  const std::vector<cv::Vec3f> &joint_pos,
				  const std::vector<float> &confidences)
{
	std::vector<int> links0 = {0, 0, 0, 0, 2, 3, 4, 5, 1, 0,  0,  9, 10, 11, 12};
	std::vector<int> links1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
	std::ofstream file (str);
	for (size_t i = 0; i < joint_pos.size (); ++i)
	{
		file << "1.0,";
		file << joint_pos[i][0] << "," << joint_pos[i][1] << "," << joint_pos[i][2] << ",";
		file << "0.0,0.0,1.0," << confidences[i] << ",";
		file << links0[i] << "," << links1[i] << ",";
		file << "0.0,0.0,0.0,1.0" << std::endl;
	}
	file.close ();
}


int 
main (int argc, 
	  char **argv)
{
	OpenNI::initialize ();

	Device mDevice;
	mDevice.open (ANY_DEVICE);

	VideoStream mDepthStream;
	mDepthStream.create (mDevice, SENSOR_DEPTH);

	VideoMode mDepthMode;
	mDepthMode.setResolution (640, 480);
	mDepthMode.setFps (30);
	mDepthMode.setPixelFormat (PIXEL_FORMAT_DEPTH_1_MM);
	mDepthStream.setVideoMode (mDepthMode);

	VideoStream mColorStream;
	mColorStream.create (mDevice, SENSOR_COLOR);
	VideoMode mColorMode;
	mColorMode.setResolution (640, 480);
	mColorMode.setFps (30);
	mColorMode.setPixelFormat (PIXEL_FORMAT_RGB888);
	mColorStream.setVideoMode (mColorMode);
	mDevice.setImageRegistrationMode (IMAGE_REGISTRATION_DEPTH_TO_COLOR);

	NiTE::initialize ();
	
	UserTracker user_tracker;
	if (user_tracker.create () != nite::STATUS_OK)
	{
		cerr << "Can't create user tracker" << endl;
		return (-1);
	}
	
	cv::namedWindow ("Depth Image", CV_WINDOW_AUTOSIZE);
	cv::namedWindow ("Color Image",  CV_WINDOW_AUTOSIZE);

	map< HandId,vector<cv::Point2f> > mapHandData;
	map< HandId,float > mapHanddepth;
	vector<cv::Point2f> vWaveList;
	vector<cv::Point2f> vClickList;
	cv::Point2f ptSize( 3, 3 );

	mDepthStream.start();
	mColorStream.start();
	

	std::vector<cv::Mat> depth_images_vec, color_images_vec, user_maps_vec;
	std::vector<std::vector<float> > confidences_vec;
	std::vector<std::vector<cv::Vec3f> > joint_pos_vec;


	int iMaxDepth = mDepthStream.getMaxPixelValue ();
	// start
	while (true)
	{
		cv::Mat cImageBGR;
		VideoFrameRef mColorFrame;
		mColorStream.readFrame (&mColorFrame);
		const cv::Mat mImageRGB (mColorFrame.getHeight (), mColorFrame.getWidth (),
			                     CV_8UC3, (void*)mColorFrame.getData ());

		// RGB ==> BGR
		cv::cvtColor (mImageRGB, cImageBGR, CV_RGB2BGR);


		UserTrackerFrameRef user_frame;
		if (user_tracker.readFrame (&user_frame) == nite::STATUS_OK)
		{
			openni::VideoFrameRef mDepthFrame = user_frame.getDepthFrame ();

			const cv::Mat mImageDepth (mDepthFrame.getHeight (), mDepthFrame.getWidth (), CV_16UC1, (void*)mDepthFrame.getData ());

			cv::Mat mScaledDepth, mImageBGR;
			mImageDepth.convertTo (mScaledDepth, CV_8U, 255.0 / iMaxDepth);

			cv::cvtColor (mScaledDepth, mImageBGR, CV_GRAY2BGR);


			const nite::Array<nite::UserData> *users = &user_frame.getUsers ();	
			int count_tracked = 0;
			int best_id = -1;
			float best_conf = 0.;
		    for (int i = 0; i < users->getSize(); ++i) 
		    {
                const nite::UserData &user = (*users)[i];
 
                if (user.isNew()) 
                    user_tracker.startSkeletonTracking(user.getId());

                if ((*users)[i].getSkeleton ().getState () == nite::SKELETON_TRACKED &&
                	(*users)[i].isVisible ())
                {
                	count_tracked ++;
                	if ((*users)[i].getSkeleton().getJoint(nite::JOINT_HEAD).getPositionConfidence() > best_conf)
					{
						best_conf = (*users)[i].getSkeleton().getJoint(nite::JOINT_HEAD).getPositionConfidence();
						best_id = i;
					}	
                }
            }
            
            /// Draw the skeleton on the depth map
            std::vector<cv::Vec3f> joint_pos;
            std::vector<float> confidences;
            if (best_id != -1)
            {
            	getJointPos ((*users)[best_id].getSkeleton (),
            				 joint_pos, confidences);

	            confidences_vec.push_back (confidences);
    	        joint_pos_vec.push_back (joint_pos);
				depth_images_vec.push_back (mImageDepth.clone ());
				color_images_vec.push_back (cImageBGR.clone ());

				visualizeSkeleton (joint_pos, cImageBGR);

				nite::UserMap user_map = user_frame.getUserMap ();
				const cv::Mat user_image (user_map.getHeight (), user_map.getWidth (), CV_16UC1, (void*)user_map.getPixels ());

				/// DEBUG
				cv::Mat user_image_manual = cv::Mat::zeros (user_image.rows, user_image.cols, CV_8UC1);
				int count_user_pixels = 0;
				for (size_t x = 0; x < user_image.cols; ++x)
					for (size_t y = 0; y < user_image.rows; ++y)
						if (user_image.at<unsigned short> (cv::Point (x, y)) != 0)
						{
							count_user_pixels ++;
							user_image_manual.at<unsigned char> (cv::Point (x, y)) = 255;
						}
				printf ("count user pixels: %d\n", count_user_pixels);

				cv::Mat aux_img;
				cv::cvtColor( user_image, aux_img, CV_GRAY2BGR );

				cv::imshow ("user map", user_image_manual);

				user_maps_vec.push_back (user_image_manual.clone ());
			}

			char str[64];
            sprintf (str, "tracked users: %d, best id %d, conf %f", count_tracked, best_id, best_conf);
            cv::putText (cImageBGR, str, cv::Point (50, 50), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar (0, 0, 255));

		
			cv::imshow( "Depth Image", mImageBGR );
			cv::imshow("Color Image", cImageBGR);

			user_frame.release();
		}
		else
			cerr << "Can't get new frame" << endl;
		
		if( cv::waitKey( 1 ) == 'q' )
			break;
	}


	printf ("saving data to hdd ...\n");
	for (size_t f_i = 0; f_i < confidences_vec.size (); ++f_i)
	{
		char str[512];
		sprintf (str, "output/%05zu_rgb.jpg", f_i);
		cv::imwrite (str, color_images_vec[f_i]);

		sprintf (str, "output/%05zu_depth.pgm", f_i);
		cv::imwrite (str, depth_images_vec[f_i]);

		sprintf (str, "output/%05zu_userMap.png", f_i);
		cv::imwrite (str, user_maps_vec[f_i]);

		sprintf (str, "output/%05zu_skel.txt", f_i);
		saveSkeletonFile (str, joint_pos_vec[f_i], confidences_vec[f_i]);
	}


	user_tracker.destroy();
	mColorStream.destroy();
	NiTE::shutdown();
	OpenNI::shutdown();

	return (0);
}
