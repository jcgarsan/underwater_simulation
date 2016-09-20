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

#ifndef OCULUSCAMERAMANIPULATOR_H
#define OCULUSCAMERAMANIPULATOR_H

#include <osgGA/CameraManipulator>
#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <osg/MatrixTransform>


class OculusCameraManipulator: public osgGA::CameraManipulator
{

	class UnknownTransformType {};

    public:
		OculusCameraManipulator(osg::MatrixTransform *node): CameraManipulator::CameraManipulator(), _mat(node)
	 	{
			ocm_sub_ = nh_.subscribe<geometry_msgs::Quaternion>("oculus/orientation", 1, &OculusCameraManipulator::ocmCallback, this);
			//Describe camera positions here;
			double *position1= new double[3];
			double *rotation1= new double[3];
			position1[0]=-0.5;position1[1]=0;position1[2]=0;
			rotation1[0]=0;rotation1[1]=90;rotation1[2]=-90;

			double *position2= new double[3];
			double *rotation2= new double[3];
			position2[0]=0.5;position2[1]=0;position2[2]=-3;
			rotation2[0]=0;rotation2[1]=180;rotation2[2]=270;

			//Add camera positions to vectors
			camPositions.push_back(position1);
			camRotations.push_back(rotation1);
			camPositions.push_back(position2);
			camRotations.push_back(rotation2);

			//select camera number
			cam=0;
		};
        virtual bool handle(const osgGA::GUIEventAdapter&, osgGA::GUIActionAdapter&);
        virtual void setByMatrix(const osg::Matrixd& m);
        virtual osg::Matrixd getMatrix() const;
        virtual void setByInverseMatrix(const osg::Matrixd &m);
        virtual osg::Matrixd getInverseMatrix() const;
		void nextCam(void);

    private:         
        osg::MatrixTransform *_mat;
        osg::Matrixd _offset;				//quaternion de oculus
		std::vector<double*> camPositions;
		std::vector<double*> camRotations;
		int cam;
		ros::NodeHandle nh_;
		ros::Subscriber ocm_sub_;
		void ocmCallback(const geometry_msgs::Quaternion::ConstPtr& hmdimu);
};

#endif
