
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

The required console arguments are the following:
* -dataset_dir -> path to the dataset captured using the *openni_nite_acquisition* - saved by default in the *output* folder
* -fbx -> path to the neutral average body with an integrated skeleton, provided in *data/skeleton.fbx*
* -bs_dir -> path to the directory containing the blendshapes for modeling the body shape. Two such sets of blendshapes are provided in *data/*
* -good_points -> path to a text file containing a list of indices for the vertices to be used in the registration and modeling process. *data/indices_* 
* -pca -> path to a pose PCA model. We provide an example in *data/motion.pca*
* -out_dir -> path to the folder where to place the output files.


**NOTE** The skeleton, neutral body mesh, and the blendshapes have been created by using the [MakeHuman Project](org). In addition to the licencing terms of the code in this repository, please respect their licensing rules too.

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