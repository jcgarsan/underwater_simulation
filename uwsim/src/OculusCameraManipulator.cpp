/* 
 * Copyright (c) 2013 University of Jaume-I.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v3.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors:
 *     Juan Carlos García
 */

#include <uwsim/OculusCameraManipulator.h>
#include <uwsim/UWSimUtils.h>


void OculusCameraManipulator::ocmCallback(const geometry_msgs::Quaternion::ConstPtr& hmdimu)
{
	//cout << "hmdimu: (" << hmdimu->x << "," << hmdimu->y << "," << hmdimu->z << "," << hmdimu->w << ")" << endl;
	_offset.makeRotate(osg::Quat(-hmdimu->x, -hmdimu->y, -hmdimu->z, hmdimu->w));
  _offset.preMultRotate(osg::Quat(osg::DegreesToRadians(camRotations[cam][2]), osg::Vec3(0,0,1)));
	_offset.preMultRotate(osg::Quat(osg::DegreesToRadians(camRotations[cam][1]), osg::Vec3(0,1,0)));
	_offset.preMultRotate(osg::Quat(osg::DegreesToRadians(camRotations[cam][0]), osg::Vec3(1,0,0)));
	_offset.preMultTranslate(osg::Vec3d(camPositions[cam][0],camPositions[cam][1],camPositions[cam][2]));
}

void OculusCameraManipulator::nextCam(void)
{
  cam++;
  cam=cam%camPositions.size();
}

osg::Matrixd OculusCameraManipulator::getMatrix() const
{
    osg::Matrix mat;          

    // get the matrix of the object node   
 	  //mat = _mat->getMatrix();
 	  mat = *getWorldCoords(_mat);

    // apply the offset = quaternion de oculus
    mat = osg::Matrixd::inverse(_offset) * mat;

    // this is where the camera should be
    return mat;
}


// we don't want the user to be able to set the matrix
void OculusCameraManipulator::setByMatrix(const osg::Matrixd& m)
{
}

void OculusCameraManipulator::setByInverseMatrix(const osg::Matrixd& m) 
{
}

// getting the inverse is easy
osg::Matrixd OculusCameraManipulator::getInverseMatrix() const 
{
    return osg::Matrix::inverse(getMatrix()); 
}

// we have nothing to handle here (no device)
bool OculusCameraManipulator::handle(const osgGA::GUIEventAdapter&, osgGA::GUIActionAdapter&) 
{
    return false; 
}


