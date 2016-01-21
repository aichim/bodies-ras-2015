#pragma once
#include "fbxsdk.h"
#include <stdio.h>

class FbxIOBasic{
protected:
    FbxManager* SDKManager;
    FbxIOSettings* ios;
    FbxScene* scene;
public:
    std::string extension(){ return ("fbx"); }
    void convertFbxdouble3ToArray(double array[3], FbxDouble3 fbx_arr){
        for(int i=0; i<3; i++){
            array[i] = (double)fbx_arr.mData[i];
        }
    }
    bool initializeFBX(const std::string& filename){
        // Create FbxManager object
        SDKManager = FbxManager::Create();

        // Create the IO settings object.
        ios = FbxIOSettings::Create(SDKManager, IOSROOT);
        SDKManager->SetIOSettings(ios);

        // Create an importer using the SDK manager.
        FbxImporter* importer = FbxImporter::Create(SDKManager,"");

        // Use the first argument as the filename for the importer.
        if(!importer->Initialize(filename.c_str(), -1, SDKManager->GetIOSettings())) {
            std::cerr << "Call to FbxImporter::Initialize() failed.";
            std::cerr << "Error returned: " << importer->GetStatus().GetErrorString();
            return (false);
        }
        else;
//            std::cerr<<"Importer Initialised";

        // Create a new scene so that it can be populated by the imported file.
        scene = FbxScene::Create(SDKManager,"myScene");

        // Import the contents of the file into the scene.
        importer->Import(scene);

        // The file is imported, so get rid of the importer.
        importer->Destroy();

        // std::cerr << "Import successful";

        return true;
    }
    bool destroyFBX(){
        // Destroy the SDK manager and all the other objects it was handling.
        SDKManager->Destroy();
        return true;
    }
};
