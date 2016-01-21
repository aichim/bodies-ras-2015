
Semantic Parametric Body Shape Estimation from Noisy Depth Sequences
====================================================================


### Dependencies (versions as tested):


* Point Cloud Library - last tested with commit `625af2db26bc6ba60ad6661f4c989bf5fe7b6a96` of [PCL repo](https://github.com/PointCloudLibrary/pcl)
* OpenCV 2.4.11
* Autodesk FBX SDK 2014 
* QGLViewer 2.6.3
* and all the other libraries necessary the compile the previous ones.



#### Acqusition - `openni_nite_acquisition`
This is a helper application used to collect RGB, depth images, as well as landmark positions from a sensor connected to the computer. It needs OpenNI and NiTE2 to work, not bundled here due to licensing issues.


#### tracking - `tracking_modeling_online`, `tracking_modeling_offline`
This is the main tracking application.


#### Generate data for rendering - `process_results`
Generate mesh sequences from the tracking results.




If you use this code, please cite the following publication:

```
@article{Ichim:2016:SPB:2873083.2873487,
 author = {Ichim, Alexandru Eugen and Tombari, Federico},
 title = {Semantic Parametric Body Shape Estimation from Noisy Depth Sequences},
 journal = {Robot. Auton. Syst.},
 issue_date = {January 2016},
 volume = {75},
 number = {PB},
 month = jan,
 year = {2016},
 issn = {0921-8890},
 pages = {539--549},
 numpages = {11},
 url = {http://dx.doi.org/10.1016/j.robot.2015.09.029},
 doi = {10.1016/j.robot.2015.09.029},
 acmid = {2873487},
 publisher = {North-Holland Publishing Co.},
 address = {Amsterdam, The Netherlands, The Netherlands},
 keywords = {3D body modeling, 3D body tracking, Depth data, Point Cloud Library},
} 
```

For any questions, remarks, corrections, please feel free to get in touch with the authors.