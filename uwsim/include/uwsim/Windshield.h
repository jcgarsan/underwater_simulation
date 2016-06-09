/* Copyright (c) 2016 University of Jaume-I.
 * 
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v3.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors:
 *      Joao Sebra
 *      Juan Carlos García
 *      Javier Perez
 */


#ifndef WINDSHIELD_H_
#define WINDSHIELD_H_

#include <iostream>
#include <string>

#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/TexMat>
#include <osg/PolygonMode>
#include <osg/NodeVisitor>
#include <osg/Node>
#include <osg/MatrixTransform>
#include <osg/PositionAttitudeTransform>
#include <osg/Point>

#include <nav_msgs/Odometry.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Int8MultiArray.h>
#include <std_msgs/Int8.h>
#include <underwater_sensor_msgs/Pressure.h>
#include <sensor_msgs/Range.h>

#include "VirtualCamera.h"
#include "ConfigXMLParser.h"
#include "VirtualRangeSensor.h"
#include "VirtualSLSProjector.h"
#include "ObjectPicker.h"
#include "InertialMeasurementUnit.h"
#include "PressureSensor.h"
#include "GPSSensor.h"
#include "DVLSensor.h"
#include "MultibeamSensor.h"
#include "Windshield.h"
#include <uwsim/SceneBuilder.h>
#include <uwsim/osgOceanScene.h>
#include <uwsim/SimulatorConfig.h>
#include <uwsim/SimulatedIAUV.h>
#include <uwsim/URDFRobot.h>

#define pressureThreshold   0.5
#define num_sensors         2       // 0 = is there an alarm?, 1 = surface, 2 = seafloor



class virtualHandleCallback : public osg::NodeCallback
{
public:
    void virtual_joy_Callback(const std_msgs::Float64MultiArray::ConstPtr& thrusterArray)
    {
        _x = thrusterArray->data[0] * 100;
        _y = thrusterArray->data[4] * 100;
        _z = 100;
        //_z = thrusterArray->data[2] * 100;
        //cout << "_x = " << _x << "; _y = " << _y << "; _z = " << _z << endl;
    }


    // Override the constructor
    virtualHandleCallback()
    {
        _x = _y = _z = 0;
        virtual_joy_sub_ = nh_.subscribe<std_msgs::Float64MultiArray>("g500/thrusters_input", 1,\
                                                 &virtualHandleCallback::virtual_joy_Callback, this);
    };
    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        float size = 0.5;
        osg::Group* currGroup = node->asGroup();
        osg::Node* foundNode;
        
        if (currGroup)
        {
            for (unsigned int i = 0 ; i < currGroup->getNumChildren(); i ++)
            {
                if (currGroup->getChild(i)->getName() == "handle_transform")
                {
                    osg::PositionAttitudeTransform* transform = dynamic_cast<osg::PositionAttitudeTransform*>(currGroup->getChild(i));
                    transform->setScale(osg::Vec3(1, 1, size));
                    
                    osg::Matrix matrix;
                    matrix.makeLookAt(osg::Vec3(0,0,0), osg::Vec3(_x,_y,_z), osg::Vec3(0,0,0));
                    
                    osg::Quat quat;
                    quat.set(matrix);
                    
                    transform->setAttitude(quat);
                }

                if (currGroup->getChild(i)->getName() == "cube_transform")
                {
                    osg::PositionAttitudeTransform* transform = dynamic_cast<osg::PositionAttitudeTransform*>(currGroup->getChild(i));
                    
                    osg::Matrix matrix;
                    matrix.makeLookAt(osg::Vec3(0,0,0), osg::Vec3(_x,_y,_z), osg::Vec3(0,0,0));
                    
                    osg::Quat quat;
                    quat.set(matrix);
                    
                    transform->setAttitude(quat);
                }
            }
        }
        
        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    float _x, _y, _z;
    float _px, _py, _pz;
    ros::NodeHandle nh_;
    ros::Subscriber virtual_joy_sub_;
};



class upDownCallback : public osg::NodeCallback
{
public:
    void robotZCallback(const nav_msgs::Odometry::ConstPtr& odom)
    {
        if(abs(_z-odom->pose.pose.position.z) < 0.01)
            _state = 0;
        else if(_z < odom->pose.pose.position.z)
            _state = 1;
        else if(_z > odom->pose.pose.position.z)
            _state = -1;    
        
        //cout << "state: " << _state << "  z: " << odom->pose.pose.position.z << "  _z: " << _z << endl;
        _z = odom->pose.pose.position.z;
    };

    // Override the constructor
    upDownCallback()
    {
        robot_z_sub_ = nh_.subscribe<nav_msgs::Odometry>("uwsim/girona500_odom_RAUVI", 1, &upDownCallback::robotZCallback, this);
    };

    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        float a,b;
        float angle;
        osg::Group* currGroup           = node->asGroup();
        osg::Geode* geode               = dynamic_cast<osg::Geode*>(currGroup->getChild(0));
        osg::Geometry* geometry         = geode->getDrawable(0)->asGeometry();
        osg::Vec4Array* colors          = new osg::Vec4Array;
        osg::MatrixTransform* transform = dynamic_cast<osg::MatrixTransform*>(node);        

        if(_state == 0)
            a = 0;
        else
            a = 1;

        if(_state == -1)
        {
            angle = 90+180;
            b = 0.7;
        }
        else if (_state == 1)
        {
            angle = 90;
            b = 1;
        }

        colors->push_back(osg::Vec4(0,0,b,a));

        geometry->setColorArray(colors);
        geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

        transform->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(angle), 0.0, 1.0, 0.0));

        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    float _z;
    int _state; //-1 -> down // 0 -> deadzone  // 1 -> up
    ros::NodeHandle nh_;
    ros::Subscriber robot_z_sub_;
};



class pressureWarningCallback : public osg::NodeCallback
{
public:
    void pressureSensorCallback(const underwater_sensor_msgs::Pressure::ConstPtr& pressureValue)
    {
        if (abs(pressureValue->pressure) < pressureThreshold)
            _state = 1;
        else
            _state = 0;
    }

    // Override the constructor
    pressureWarningCallback()
    {
        trans_state   = 0;
        _transparency = 1;
        _state        = 0;
        snsrPressure_sub_ = nh_.subscribe<underwater_sensor_msgs::Pressure>("g500/pressure", 1, &pressureWarningCallback::pressureSensorCallback, this);
    };

    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        float alphaChannel;
        osg::Group* currGroup   = node->asGroup();
        osg::Geode* geode       = dynamic_cast<osg::Geode*>(currGroup->getChild(0));
        osg::Geometry* geometry = geode->getDrawable(0)->asGeometry();
        osg::Vec4Array* colors  = new osg::Vec4Array;

        if(trans_state)
            _transparency += 0.2;
        else
            _transparency -= 0.2;

        if(_transparency >= 1)
            trans_state = 0;
        else if(_transparency <=0)
            trans_state = 1;

        if(_state == 0)
            alphaChannel = 0;
        else
            alphaChannel = _transparency;

        colors->push_back(osg::Vec4(1, 1, 1, alphaChannel));

        geometry->setColorArray(colors);
        geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    float           _transparency;
    bool            trans_state;
    int             _state; //-1 -> down // 0 -> deadzone  // 1 -> up
    ros::NodeHandle nh_;
    ros::Subscriber snsrPressure_sub_;
};



class warningSafetyAlarmCallback : public osg::NodeCallback
{
public:
    void warningCallback(const std_msgs::Int8MultiArray::ConstPtr& msg)
    {
        _state = msg->data[0];
    }

    // Override the constructor
    warningSafetyAlarmCallback()
    {
        trans_state   = 0;
        _transparency = 1;
        _state        = 0;
        warning_sub_ = nh_.subscribe<std_msgs::Int8MultiArray>("safetyMeasuresAlarm", 1, &warningSafetyAlarmCallback::warningCallback, this);
    };

    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        float alphaChannel;
        osg::Group* currGroup   = node->asGroup();
        osg::Geode* geode       = dynamic_cast<osg::Geode*>(currGroup->getChild(0));
        osg::Geometry* geometry = geode->getDrawable(0)->asGeometry();
        osg::Vec4Array* colors  = new osg::Vec4Array;

        if(trans_state)
            _transparency += 0.2;
        else
            _transparency -= 0.2;

        if(_transparency >= 1)
            trans_state = 0;
        else if(_transparency <=0)
            trans_state = 1;

        if(_state == 0)
            alphaChannel = 0;
        else
            alphaChannel = _transparency;

        colors->push_back(osg::Vec4(1, 1, 1, alphaChannel));

        geometry->setColorArray(colors);
        geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    float               _transparency;
    bool                trans_state;
    int                 _state; //-1 -> down // 0 -> deadzone  // 1 -> up
    ros::NodeHandle     nh_;
    ros::Subscriber     warning_sub_;
};



class userControlCallback : public osg::NodeCallback
{
public:
    void userCallback(const std_msgs::Int8MultiArray::ConstPtr& msg)
    {
        _state = msg->data[0];
    }

    // Override the constructor
    userControlCallback()
    {
        trans_state   = 0;
        _transparency = 1;
        _state        = 0;
        userControl_sub_ = nh_.subscribe<std_msgs::Int8MultiArray>("userControlAlarm", 1, &userControlCallback::userCallback, this);
    };

    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        float alphaChannel;
        osg::Group* currGroup   = node->asGroup();
        osg::Geode* geode       = dynamic_cast<osg::Geode*>(currGroup->getChild(0));
        osg::Geometry* geometry = geode->getDrawable(0)->asGeometry();
        osg::Vec4Array* colors  = new osg::Vec4Array;

        if(trans_state)
            _transparency += 0.2;
        else
            _transparency -= 0.2;

        if(_transparency >= 1)
            trans_state = 0;
        else if(_transparency <=0)
            trans_state = 1;

        if(_state == 0)
            alphaChannel = 0;
        else
            alphaChannel = _transparency;

        colors->push_back(osg::Vec4(1, 1, 1, alphaChannel));

        geometry->setColorArray(colors);
        geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    float               _transparency;
    bool                trans_state;
    int                 _state; //-1 -> down // 0 -> deadzone  // 1 -> up
    ros::NodeHandle     nh_;
    ros::Subscriber     userControl_sub_;
};



class missionControlCallback : public osg::NodeCallback
{
public:
    void missionCallback(const std_msgs::Int8::ConstPtr& msg)
    {
        _state = msg->data;
    }

    // Override the constructor
    missionControlCallback()
    {
        _state = -1;
        missionControl_sub_ = nh_.subscribe<std_msgs::Int8>("missionControlAlarm", 1, &missionControlCallback::missionCallback, this);
    };


    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        float alphaChannel;
        osg::Group* currGroup   = node->asGroup();
        osg::Geode* geode       = dynamic_cast<osg::Geode*>(currGroup->getChild(0));
        osgText::Text* text     = (osgText::Text *) geode->getDrawable(0);;

        if (_state == 1)
        {
            alphaChannel = 1;
            text->setText("Mission succeed");
        }
        else if (_state == 0)
        {
            alphaChannel = 1;
            text->setText("Mission failed");
        }
        else
        {
            alphaChannel = 0;
            text->setText("");
        }            

        text->setColor(osg::Vec4(1, 1, 1, alphaChannel));

        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    int                 _state; //-1 -> mission in progress // 0 -> fail  // 1 -> success
    ros::NodeHandle     nh_;
    ros::Subscriber     missionControl_sub_;
};



class RotationColorCallback : public osg::NodeCallback
{
public:
    // Override the constructor
    RotationColorCallback()
    {
        _angle = -145;
        state = 0;
    };


    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        osg::MatrixTransform* transform = dynamic_cast<osg::MatrixTransform*>(node);

        if (!state)
        {
            if (_angle > -145)
                _angle -= 1;
            else
                state = 1;
        }
        else
        {
            if (_angle < 150)
                _angle += 1;
            else
                state = 0;
        }
        
        transform->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(_angle), 0.0, 1.0, 0.0));
        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    bool state;
    float _angle;
};



class sliderGaugeCallback : public osg::NodeCallback
{
public:
    void sensorRngCallback(const sensor_msgs::Range::ConstPtr& rangeValue)
    {
        _percentage = 1-(rangeValue->range) / (rangeValue->max_range);
    }

    // Override the constructor
    sliderGaugeCallback(float width, float height)
    {
        _width          = width;
        _height         = height;
        _percentage     = 0;
        state           = 0;
        snsrRange_sub_  = nh_.subscribe<sensor_msgs::Range>("uwsim/g500/range", 1, &sliderGaugeCallback::sensorRngCallback, this);
    };
    
    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        osg::Group* currGroup = node->asGroup();
        osg::Node* foundNode;
        float r, g, a;
        
        if (currGroup) 
        {
            for (unsigned int i = 0 ; i < currGroup->getNumChildren(); i ++)
            {
                if (currGroup->getChild(i)->getName() == "bar_status")
                {
                    osg::Geode* geode = dynamic_cast<osg::Geode*>(currGroup->getChild(i));
                    osg::Geometry* geometry = geode->getDrawable(0)->asGeometry();

                    r = _percentage;
                    g = 1.0 - _percentage;
                    a = 1.2 - _percentage;
                    
                    osg::Vec4Array* colors = new osg::Vec4Array;
                    colors->push_back(osg::Vec4(r,g,0,a));

                    geometry->setColorArray(colors);
                    geometry->setColorBinding(osg::Geometry::BIND_OVERALL);
                    
                    osg::Vec3Array* vertices = new osg::Vec3Array(5);
                    (*vertices)[0].set(-_width/2, 0, -_height/2);
                    (*vertices)[1].set(-_width/2, 0, _height/2 - (1-_percentage)*_height);
                    (*vertices)[2].set(_width/2, 0, _height/2 - (1-_percentage)*_height);
                    (*vertices)[3].set(_width/2, 0, -_height/2);
                    (*vertices)[4].set(-_width/2, 0, -_height/2);
                    
                    // pass the created vertex array to the points geometry object.
                    geometry->setVertexArray(vertices);
                }
            }
        }
        
        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };

private:
    bool            state;
    float           _percentage;
    float           _width, _height;
    ros::NodeHandle nh_;
    ros::Subscriber snsrRange_sub_;
};



osg::Geode* createSphere(float radius, osg::Vec4 color)
{
    osg::Sphere* unitCube                   = new osg::Sphere(osg::Vec3(0,0,0), radius);
    osg::ShapeDrawable* unitCubeDrawable    = new osg::ShapeDrawable(unitCube);
    osg::Geode* basicShapesGeode            = new osg::Geode();
    osg::ref_ptr<osg::BlendFunc> blendFunc  = new osg::BlendFunc;

    unitCubeDrawable->setColor(color);

    basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON); //Unset shader

    blendFunc->setFunction(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(blendFunc);

    basicShapesGeode->addDrawable(unitCubeDrawable);
    return basicShapesGeode;
}



osg::Group* createVirtualHandle()
{
    osg::Cylinder* cylinder = new osg::Cylinder(osg::Vec3(0,0,0.5),0.05,1);
    osg::ShapeDrawable* cylinderDrawable = new osg::ShapeDrawable(cylinder);
    osg::Geode* cylinderGeode = new osg::Geode();
    cylinderGeode->addDrawable(cylinderDrawable);
    cylinderGeode->setName("cylinder");

    osg::PositionAttitudeTransform* transform = new osg::PositionAttitudeTransform;
    transform->addChild(cylinderGeode);
    transform->setName("handle_transform");
    
    osg::Geode* cube = createSphere(0.15, osg::Vec4(0.8,0.8,0.8,1));

    osg::PositionAttitudeTransform* transform2 = new osg::PositionAttitudeTransform;
    transform2->addChild(cube);
    transform2->setName("cube_transform");
    
    osg::Group* group = new osg::Group;
    
    group->addChild(transform);
    group->addChild(transform2);
    group->setUpdateCallback(new virtualHandleCallback());
    
    return group;
}



osg::Geode* drawCube(int cx, int cy, int cz, float size)
{
    osg::Box* unitCube = new osg::Box( osg::Vec3(cx,cy,cz), size);
    osg::ShapeDrawable* unitCubeDrawable = new osg::ShapeDrawable(unitCube);
    osg::Geode* basicShapesGeode = new osg::Geode();
    basicShapesGeode->addDrawable(unitCubeDrawable);

    return basicShapesGeode;
}



osg::Geode* createBox(osg::Vec3 center, float tx, float ty, float tz, osg::Vec4 color)
{
    osg::Box* unitCube = new osg::Box(center, -tx, -ty, -tz);
    osg::ShapeDrawable* unitCubeDrawable = new osg::ShapeDrawable(unitCube);

    unitCubeDrawable->setColor(color);

    osg::Geode* basicShapesGeode = new osg::Geode();

    basicShapesGeode->addDrawable(unitCubeDrawable);
    return basicShapesGeode;
}



osg::Geode* createQuadWithTex(float size, float r, float g, float b, float a, const std::string& filename)
{
    osg::ref_ptr<osg::Vec4Array> shared_colors = new osg::Vec4Array;
    shared_colors->push_back(osg::Vec4(r,g,b,a));

    osg::ref_ptr<osg::Vec3Array> shared_normals = new osg::Vec3Array;
    shared_normals->push_back(osg::Vec3(0.0f,-1.0f,0.0f));

    // create Geometry object to store all the vertices and lines primitive.
    osg::Geometry* polyGeom = new osg::Geometry();

    // note, anticlockwise ordering.
    osg::Vec3 myCoords[] =
    {
        osg::Vec3(-size/2, 0, -size/2),
        osg::Vec3(-size/2, 0, size/2),
        osg::Vec3(size/2, 0, size/2),
        osg::Vec3(size/2, 0, -size/2),
    };

    int numCoords = sizeof(myCoords)/sizeof(osg::Vec3);

    osg::Vec3Array* vertices = new osg::Vec3Array(numCoords,myCoords);

    // pass the created vertex array to the points geometry object.
    polyGeom->setVertexArray(vertices);
   
    // use the shared color array.
    polyGeom->setColorArray(shared_colors);
    polyGeom->setColorBinding(osg::Geometry::BIND_OVERALL);

    // use the shared normal array.
    polyGeom->setNormalArray(shared_normals);
    polyGeom->setNormalBinding(osg::Geometry::BIND_OVERALL);


    // This time we simply use primitive, and hardwire the number of coords to use
    // since we know up front,
    polyGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS,0,numCoords));

    osg::Vec2Array* texcoords = new osg::Vec2Array(4);
    (*texcoords)[0].set(0.0f, 0.0f);
    (*texcoords)[1].set(1.0f, 0.0f);
    (*texcoords)[2].set(1.0f, 1.0f);
    (*texcoords)[3].set(0.0f, 1.0f);
    polyGeom->setTexCoordArray(0,texcoords);
    polyGeom->setUseDisplayList(false);

    // load image
    osg::Image* img = osgDB::readImageFile(filename);

    // setup texture
    osg::TextureRectangle* texture = new osg::TextureRectangle(img);
    osg::TexMat* texmat = new osg::TexMat;
    texmat->setScaleByTextureRectangleSize(true);

    osg::Geode* polygeode = new osg::Geode();

    polygeode->getOrCreateStateSet()->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);

    // setup state
    osg::StateSet* state = polyGeom->getOrCreateStateSet();
    state->setTextureAttributeAndModes(0, texture, osg::StateAttribute::ON);
    state->setTextureAttributeAndModes(0, texmat, osg::StateAttribute::ON);

    polygeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON); //Unset shader

    // turn off lighting
    state->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

    osg::ref_ptr<osg::BlendFunc> blendFunc = new osg::BlendFunc;
    blendFunc->setFunction(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    polygeode->getOrCreateStateSet()->setAttributeAndModes(blendFunc);
    polygeode->addDrawable(polyGeom);
    polygeode->setName("texture_geode");

    return polygeode;
}



osg::Geode* createTextBox(float radius, const std::string& text)
{
    osg::Geode* geode  = new osg::Geode;
    float characterSize = radius*0.2f;

    osgText::Text* text1 = new osgText::Text;
    text1->setFont("~/.uwsim/data/objects/arial.ttf");
    text1->setCharacterSize(characterSize);
    text1->setPosition(osg::Vec3(0,0,0));
    text1->setAxisAlignment(osgText::Text::XZ_PLANE);
    text1->setColor(osg::Vec4(1,0,0,1));
    text1->setText(text);
    text1->setAlignment(osgText::Text::CENTER_CENTER);
    text1->setDrawMode(osgText::Text::TEXT|osgText::Text::ALIGNMENT|osgText::Text::BOUNDINGBOX);
    text1->setDrawMode(osgText::Text::TEXT);
    //text1->setDataVariance(DYNAMIC); 
	text1->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON);
    geode->addDrawable(text1);

    return geode;
}



osg::Group* create3DText(const osg::Vec3& center, const std::string& originalText)
{

    osg::Group* rootNode = new osg::Group;
    osg::Geode* geode    = new osg::Geode;
    osgText::Text* text  = new osgText::Text;
    //float characterSize = 0.2f;

    text->setFont("~/.uwsim/data/objects/arial.ttf");
    text->setCharacterSize(0.2f);
    text->setPosition(center);
    text->setAxisAlignment(osgText::Text::SCREEN);
    text->setText(originalText);
    geode->addDrawable(text);
    rootNode->addChild(geode);

    return rootNode;    
}



osg::Group* create3DText2()
{

    osg::Group* rootNode = new osg::Group;
    osg::Geode* geode    = new osg::Geode;
    osgText::Text* text  = new osgText::Text;
    float characterSize = 0.2f;

    text->setFont("~/.uwsim/data/objects/arial.ttf");
    text->setCharacterSize(characterSize);
    text->setPosition(osg::Vec3(0.7, -0.5, 0));
    text->setAxisAlignment(osgText::Text::SCREEN);
    text->setText("Mission succeed");
    geode->addDrawable(text);
    rootNode->addChild(geode);

    return rootNode;
}



osg::MatrixTransform* setNodePosition(float up, float right, float depth, osg::Node* geode)
{
    osg::MatrixTransform* transform = new osg::MatrixTransform;
    transform->setMatrix(osg::Matrix::translate(0, 1+depth, 0));
    transform->addChild(geode);

    osg::MatrixTransform* invert = new osg::MatrixTransform;
    osg::Vec3d axis(0, 1, 0);
    invert->setMatrix(osg::Matrix::rotate(osg::PI, axis));
    invert->addChild(transform);

    osg::MatrixTransform* transform2 = new osg::MatrixTransform;
    axis.set(1, 0, 0);
    transform2->setMatrix(osg::Matrix::rotate(up, axis));
    transform2->addChild(invert);

    osg::MatrixTransform* transform3 = new osg::MatrixTransform;
    axis.set(0, 0, 1);
    transform3->setMatrix(osg::Matrix::rotate(right, axis));
    transform3->addChild(transform2);

    return transform3;
}



osg::Geode* createRectangularBox(float width, float height, float line_width, float r, float g, float b, float a)
{
    osg::Geode* geode = new osg::Geode;
    
    // create Geometry object to store all the vertices and lines primitive.
    osg::Geometry* linesGeom = new osg::Geometry();
    
    osg::Vec3Array* vertices = new osg::Vec3Array(5);
    (*vertices)[0].set(-width/2, 0, -height/2);
    (*vertices)[1].set(-width/2, 0, height/2);
    (*vertices)[2].set(width/2, 0, height/2);
    (*vertices)[3].set(width/2, 0, -height/2);
    (*vertices)[4].set(-width/2, 0, -height/2);
    
    // pass the created vertex array to the points geometry object.
    linesGeom->setVertexArray(vertices);
    
    // set the colors as before, plus using the above
    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(r,g,b,a));
    linesGeom->setColorArray(colors);
    linesGeom->setColorBinding(osg::Geometry::BIND_OVERALL);

  // set the normal in the same way color.
    osg::Vec3Array* normals = new osg::Vec3Array;
    normals->push_back(osg::Vec3(0.0f,-1.0f,0.0f));
    linesGeom->setNormalArray(normals);
    linesGeom->setNormalBinding(osg::Geometry::BIND_OVERALL);
    
    // This time we simply use primitive, and hardwire the number of coords to use
    // since we know up front,
    linesGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_STRIP,0,5));
    
    osg::LineWidth* linewidth = new osg::LineWidth();
    linewidth->setWidth(line_width);
    linesGeom->getOrCreateStateSet()->setAttribute(linewidth, osg::StateAttribute::ON);
    
    // add the points geometry to the geode.
    geode->addDrawable(linesGeom);
    
    return geode;
}



osg::Group* createAnalogDial(float up, float right, float size, const std::string& pointer_filename, const std::string& background_filename)
{
    osg::Group* dial = new osg::Group;

    osg::Geode* background_geode = createQuadWithTex(size, 1, 1, 1, 0.5, background_filename);
    osg::Geode* pointer_geode = createQuadWithTex(size, 0, 0, 0, 0.5, pointer_filename);
    
    background_geode->setName("background_geode");
    pointer_geode->setName("pointer_geode");
    
    osg::MatrixTransform* background = setNodePosition(up, right, 0, background_geode);
    osg::MatrixTransform* pointer = new osg::MatrixTransform;

    pointer->addChild(pointer_geode);
    pointer->setUpdateCallback(new RotationColorCallback());
    
    dial->addChild(background);
    dial->addChild(setNodePosition(up, right, -0.01, pointer));
    dial->addChild(setNodePosition(up, right, 0, createRectangularBox(size+0.1, size+0.1, 1, 0, 0, 0, 0.5)));
    //dial->setUpdateCallback(new RotationColorCallback());
    
    return dial;
}



osg::Group* createBarGauge(float width, float height, float percentage)
{
    float diff = height/10;
    
    width -= diff;
    height -= diff;
    
    osg::Geode* geode = new osg::Geode;
    
    // create Geometry object to store all the vertices and lines primitive.
    osg::Geometry* polyGeom = new osg::Geometry();
    
    osg::Vec3Array* vertices = new osg::Vec3Array(5);
    (*vertices)[0].set(-width/2, 0, -height/2);
    (*vertices)[1].set(-width/2, 0, height/2 - (1-percentage)*height);
    (*vertices)[2].set(width/2, 0, height/2 - (1-percentage)*height);
    (*vertices)[3].set(width/2, 0, -height/2);
    (*vertices)[4].set(-width/2, 0, -height/2);
    
    // pass the created vertex array to the points geometry object.
    polyGeom->setVertexArray(vertices);
    
    // set the colors as before, plus using the above
    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1,1,1,1));
    polyGeom->setColorArray(colors);
    polyGeom->setColorBinding(osg::Geometry::BIND_OVERALL);
    
    // set the normal in the same way color.
    osg::Vec3Array* normals = new osg::Vec3Array;
    normals->push_back(osg::Vec3(0.0f,-1.0f,0.0f));
    polyGeom->setNormalArray(normals);
    polyGeom->setNormalBinding(osg::Geometry::BIND_OVERALL);
    
    // This time we simply use primitive, and hardwire the number of coords to use
    // since we know up front,
    polyGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS,0,5));
    
    // setup state
    osg::StateSet* state = polyGeom->getOrCreateStateSet();
    
    // turn off lighting
    state->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
    
    osg::ref_ptr<osg::BlendFunc> blendFunc = new osg::BlendFunc;
    blendFunc->setFunction(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    polyGeom->getOrCreateStateSet()->setAttributeAndModes(blendFunc);

    polyGeom->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON); //Unset shader
    
    // add the points geometry to the geode.
    geode->addDrawable(polyGeom);
    
    geode->setName("bar_status");
    
    width += diff;
    height += diff;
    
    osg::Group* group = new osg::Group;
    
    group->addChild(geode);
    group->addChild(createRectangularBox(width, height, 0, 0, 0, 1, 1));
    
    return group;
}



osg::Geode* createHUDButton(int buttonID, const osg::Vec3& center, const std::string& textButton)
{
    // A geometry node for our HUD:
    osg::Geode* HUDGeode2 = new osg::Geode();
    // Setup the button ID
    std::string Id = "button" + buttonID;
    HUDGeode2->setName(Id);
    // Text instance that wil show up in the HUDButton:
    osgText::Text* textOne = new osgText::Text();

    // Set up geometry for the HUD and add it to the HUD
    osg::Geometry* HUDBackgroundGeometry2 = new osg::Geometry();

    // Set up the button size & position
    osg::Vec3Array* HUDBackgroundVertices2 = new osg::Vec3Array;
    if (buttonID == 0)
    {   
        HUDBackgroundVertices2->push_back(osg::Vec3(678,   0, -1));
        HUDBackgroundVertices2->push_back(osg::Vec3(791,   0, -1));
        HUDBackgroundVertices2->push_back(osg::Vec3(791, 100, -1));
        HUDBackgroundVertices2->push_back(osg::Vec3(678, 100, -1));
    }
    else
    {
        HUDBackgroundVertices2->push_back(osg::Vec3(113*(buttonID-1),       0, -1));
        HUDBackgroundVertices2->push_back(osg::Vec3(113*(buttonID-1)+100,   0, -1));
        HUDBackgroundVertices2->push_back(osg::Vec3(113*(buttonID-1)+100, 100, -1));
        HUDBackgroundVertices2->push_back(osg::Vec3(113*(buttonID-1),     100, -1));
    }

    osg::DrawElementsUInt* HUDBackgroundIndices2 = new osg::DrawElementsUInt(osg::PrimitiveSet::POLYGON, 0);
    HUDBackgroundIndices2->push_back(0);
    HUDBackgroundIndices2->push_back(1);
    HUDBackgroundIndices2->push_back(2);
    HUDBackgroundIndices2->push_back(3);

    osg::Vec4Array* HUDcolors2 = new osg::Vec4Array;
    HUDcolors2->push_back(osg::Vec4(0.8f,0.8f,0.8f,0.8f));

    osg::Vec2Array* texcoords2 = new osg::Vec2Array(4);
    (*texcoords2)[0].set(0.0f,0.0f);
    (*texcoords2)[1].set(1.0f,0.0f);
    (*texcoords2)[2].set(1.0f,1.0f);
    (*texcoords2)[3].set(0.0f,1.0f);

    HUDBackgroundGeometry2->setTexCoordArray(0,texcoords2);
    osg::Texture2D* HUDTexture2 = new osg::Texture2D;
    HUDTexture2->setDataVariance(osg::Object::DYNAMIC);
    osg::Image* hudImage2;
    hudImage2 = osgDB::readImageFile("~/.uwsim/data/textures/HUD_button_background.jpg");
    HUDTexture2->setImage(hudImage2);
    osg::Vec3Array* HUDnormals2 = new osg::Vec3Array;
    HUDnormals2->push_back(osg::Vec3(0.0f,0.0f,1.0f));
    HUDBackgroundGeometry2->setNormalArray(HUDnormals2);
    HUDBackgroundGeometry2->setNormalBinding(osg::Geometry::BIND_OVERALL);
    HUDBackgroundGeometry2->addPrimitiveSet(HUDBackgroundIndices2);
    HUDBackgroundGeometry2->setVertexArray(HUDBackgroundVertices2);
    HUDBackgroundGeometry2->setColorArray(HUDcolors2);
    HUDBackgroundGeometry2->setColorBinding(osg::Geometry::BIND_OVERALL);

    HUDGeode2->addDrawable(HUDBackgroundGeometry2);

     // Create and set up a state set using the texture from above:
    osg::StateSet* HUDStateSet2 = new osg::StateSet();
    HUDGeode2->setStateSet(HUDStateSet2);
    HUDStateSet2->setTextureAttributeAndModes(0,HUDTexture2,osg::StateAttribute::ON);

    // For this state set, turn blending on (so alpha texture looks right)
    HUDStateSet2->setMode(GL_BLEND,osg::StateAttribute::ON);

    // Disable depth testing so geometry is draw regardless of depth values
    // of geometry already draw.
    HUDStateSet2->setMode(GL_DEPTH_TEST,osg::StateAttribute::OFF);
    HUDStateSet2->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);

    // Need to make sure this geometry is draw last. RenderBins are handled
    // in numerical order so set bin number to 11
    HUDStateSet2->setRenderBinDetails(11, "RenderBin");


    // Add the text (Text class is derived from drawable) to the geode:
    HUDGeode2->addDrawable(textOne);

    // Set up the parameters for the text we'll add to the HUD:
    textOne->setCharacterSize(20);
    textOne->setFont("~/.uwsim/data/objects/arial.ttf");
    textOne->setText(textButton);
    textOne->setAxisAlignment(osgText::Text::SCREEN);
//    textOne->setPosition(osg::Vec3(30, 40, 0));
    textOne->setPosition(center);
    textOne->setColor(osg::Vec4(0, 0, 0, 1));


    return HUDGeode2;
}

osg::Group* createHUD()
{
    // Initialize root of scene:
    osg::Group* root = new osg::Group();;
    // A geometry node for our HUD:
    osg::Geode* HUDGeode = new osg::Geode();
    // Text instance that wil show up in the HUD:
    osgText::Text* menuTitle = new osgText::Text();
    // Projection node for defining view frustrum for HUD:
    osg::Projection* HUDProjectionMatrix = new osg::Projection;

    // Initialize the projection matrix for viewing everything we
    // will add as descendants of this node. Use screen coordinates
    // to define the horizontal and vertical extent of the projection
    // matrix. Positions described under this node will equate to
    // pixel coordinates.
    HUDProjectionMatrix->setMatrix(osg::Matrix::ortho2D(0,1280,0,800));

    // For the HUD model view matrix use an identity matrix:
    osg::MatrixTransform* HUDModelViewMatrix = new osg::MatrixTransform;
    HUDModelViewMatrix->setMatrix(osg::Matrix::identity());

    // Make sure the model view matrix is not affected by any transforms
    // above it in the scene graph:
    HUDModelViewMatrix->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

    // Add the HUD projection matrix as a child of the root node
    // and the HUD model view matrix as a child of the projection matrix
    // Anything under this node will be viewed using this projection matrix
    // and positioned with this model view matrix.
    root->addChild(HUDProjectionMatrix);
    HUDProjectionMatrix->addChild(HUDModelViewMatrix);

    // Add the Geometry node to contain HUD geometry as a child of the
    // HUD model view matrix.
    HUDModelViewMatrix->addChild(HUDGeode);

    // Set up geometry for the HUD and add it to the HUD
    osg::Geometry* HUDBackgroundGeometry = new osg::Geometry();

    // ToDo:
    // We need to increase the Y-value above 100 pixels
    osg::Vec3Array* HUDBackgroundVertices = new osg::Vec3Array;
    HUDBackgroundVertices->push_back(osg::Vec3( 0,     0, -1));
    HUDBackgroundVertices->push_back(osg::Vec3(1280,   0, -1));
    HUDBackgroundVertices->push_back(osg::Vec3(1280, 150, -1));
    HUDBackgroundVertices->push_back(osg::Vec3(   0, 150, -1));

    osg::DrawElementsUInt* HUDBackgroundIndices = new osg::DrawElementsUInt(osg::PrimitiveSet::POLYGON, 0);
    HUDBackgroundIndices->push_back(0);
    HUDBackgroundIndices->push_back(1);
    HUDBackgroundIndices->push_back(2);
    HUDBackgroundIndices->push_back(3);

    osg::Vec4Array* HUDcolors = new osg::Vec4Array;
    HUDcolors->push_back(osg::Vec4(0.8f, 0.8f, 0.8f, 0.8f));

    osg::Vec2Array* texcoords = new osg::Vec2Array(4);
    (*texcoords)[0].set(0.0f, 0.0f);
    (*texcoords)[1].set(1.0f, 0.0f);
    (*texcoords)[2].set(1.0f, 1.0f);
    (*texcoords)[3].set(0.0f, 1.0f);

    HUDBackgroundGeometry->setTexCoordArray(0,texcoords);
    osg::Texture2D* HUDTexture = new osg::Texture2D;
    HUDTexture->setDataVariance(osg::Object::DYNAMIC);
    osg::Image* hudImage;
    hudImage = osgDB::readImageFile("~/.uwsim/data/textures/HUD_background.jpg");
    HUDTexture->setImage(hudImage);
    osg::Vec3Array* HUDnormals = new osg::Vec3Array;
    HUDnormals->push_back(osg::Vec3(0.0f, 0.0f, 1.0f));
    HUDBackgroundGeometry->setNormalArray(HUDnormals);
    HUDBackgroundGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
    HUDBackgroundGeometry->addPrimitiveSet(HUDBackgroundIndices);
    HUDBackgroundGeometry->setVertexArray(HUDBackgroundVertices);
    HUDBackgroundGeometry->setColorArray(HUDcolors);
    HUDBackgroundGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);

    HUDGeode->addDrawable(HUDBackgroundGeometry);

    // Create and set up a state set using the texture from above:
    osg::StateSet* HUDStateSet = new osg::StateSet();
    HUDGeode->setStateSet(HUDStateSet);
    HUDStateSet->setTextureAttributeAndModes(0, HUDTexture, osg::StateAttribute::ON);

    // For this state set, turn blending on (so alpha texture looks right)
    HUDStateSet->setMode(GL_BLEND, osg::StateAttribute::ON);

    // Disable depth testing so geometry is draw regardless of depth values
    // of geometry already draw.
    HUDStateSet->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF);
    HUDStateSet->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);

    // Need to make sure this geometry is draw last. RenderBins are handled
    // in numerical order so set bin number to 11
    HUDStateSet->setRenderBinDetails(11, "RenderBin");

    // Add the text (Text class is derived from drawable) to the geode:
    HUDGeode->addDrawable(menuTitle);

    // Set up the parameters for the text we'll add to the HUD:
    menuTitle->setCharacterSize(25);
    menuTitle->setFont("~/.uwsim/data/objects/arial.ttf");
    menuTitle->setText("New user HUDMenu");
    menuTitle->setAxisAlignment(osgText::Text::SCREEN);
    menuTitle->setPosition(osg::Vec3(30, 110, 0));
    menuTitle->setColor(osg::Vec4(0, 0, 0, 1));

    HUDModelViewMatrix->addChild(createHUDButton(1, osg::Vec3( 15, 40, 0), "Survey"));
    HUDModelViewMatrix->addChild(createHUDButton(2, osg::Vec3(113, 50, 0), "Object\nRecovery"));
    HUDModelViewMatrix->addChild(createHUDButton(3, osg::Vec3(226, 50, 0), "Panel\nIntervention"));
    HUDModelViewMatrix->addChild(createHUDButton(4, osg::Vec3(339, 50, 0), "Dredging\nIntervention"));
    HUDModelViewMatrix->addChild(createHUDButton(5, osg::Vec3(452, 50, 0), "Go to\nsurface"));
    HUDModelViewMatrix->addChild(createHUDButton(6, osg::Vec3(565, 50, 0), "Test the\nsystem"));
    HUDModelViewMatrix->addChild(createHUDButton(0, osg::Vec3(678, 40, 0), "Exit"));

    return root;
}




void createWindshield (osg::MatrixTransform *baseTransform)
{
    osg::PositionAttitudeTransform* transform_ = new osg::PositionAttitudeTransform;
    osg::Group* interface = new osg::Group;

    transform_->setPosition(osg::Vec3(0.4, 0, 0.7));
    transform_->addChild(interface);

    osg::Group* bar = createBarGauge(0.05, 0.4, 0.2);
    bar->setUpdateCallback(new sliderGaugeCallback(0.05, 0.4));
    transform_->addChild(setNodePosition(0.2, -osg::PI/2-0.4, 0, bar));

    osg::PositionAttitudeTransform* transform_vh = new osg::PositionAttitudeTransform;
    transform_vh->setPosition(osg::Vec3(0.7, 0.5, 0.6));    
    transform_vh->addChild(createVirtualHandle());
    transform_->addChild(transform_vh);

/*    for (float b=0; b>-1*osg::PI/2; b-=10)
    {
        transform_->addChild(setNodePosition(b, -osg::PI/2+0.5, 0, createTextBox(1 , "windshield")));
    }
*/

    osg::Group* missionStatusText = create3DText(osg::Vec3(0.7, -0.5, 0), "");
    missionStatusText->setUpdateCallback(new missionControlCallback());
    transform_->addChild(setNodePosition(0.7, -0.5, 0, missionStatusText));

    osg::MatrixTransform* udArrow =  new osg::MatrixTransform;
    udArrow->addChild(createQuadWithTex(0.2,1,1,1,1,"~/.uwsim/data/objects/up.png"));
    udArrow->setUpdateCallback(new upDownCallback());
    osg::MatrixTransform* udArrowTransform = setNodePosition(-0.1, -osg::PI/2+0.4,0,udArrow);
    transform_->addChild(udArrowTransform);

/*    osg::MatrixTransform* pressureWarning =  new osg::MatrixTransform;
    pressureWarning->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(-90.0), 0.0, 1.0, 0.0));
    pressureWarning->addChild(createQuadWithTex(0.2,1,1,1,1,"~/.uwsim/data/objects/warning.png"));
    pressureWarning->setUpdateCallback(new pressureWarningCallback());
    osg::MatrixTransform* pressureWarningTransform = setNodePosition(-0.2, -osg::PI/2,0,pressureWarning);
    transform_->addChild(pressureWarningTransform);
*/

    osg::MatrixTransform* safetyWarning =  new osg::MatrixTransform;
    safetyWarning->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(-90.0), 0.0, 1.0, 0.0));
    safetyWarning->addChild(createQuadWithTex(0.2,1,1,1,1,"~/.uwsim/data/objects/warning.png"));
    safetyWarning->setUpdateCallback(new warningSafetyAlarmCallback());
    osg::MatrixTransform* safetyWarningTransform = setNodePosition(-0.2, -osg::PI/2, 0, safetyWarning);
    transform_->addChild(safetyWarningTransform);


    osg::MatrixTransform* userControlWarning =  new osg::MatrixTransform;
    userControlWarning->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(-90.0), 0.0, 1.0, 0.0));
    userControlWarning->addChild(createQuadWithTex(0.2,1,1,1,1,"~/.uwsim/data/objects/user.png"));
    userControlWarning->setUpdateCallback(new userControlCallback());
    osg::MatrixTransform* userControlWarningTransform = setNodePosition(-0.2, -osg::PI/2+0.2, 0, userControlWarning);
    transform_->addChild(userControlWarningTransform);

    baseTransform->addChild(createBox(osg::Vec3(1.1,0,1.3),0.5,1.5,0.07,osg::Vec4(0.2,0.2,0.2,1)));
    baseTransform->addChild(transform_);

    osg::Group* menuHUD = createHUD();
    transform_->addChild(menuHUD);
}

#endif
