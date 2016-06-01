#ifndef WINDSHIELD_H_
#define WINDSHIELD_H_



#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/TexMat>
#include <osg/PolygonMode>
#include <osg/NodeVisitor>
#include <osg/Node>
#include <nav_msgs/Odometry.h>
#include <std_msgs/Float64MultiArray.h>
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
#include <osg/MatrixTransform>
#include <osg/PositionAttitudeTransform>

#include <osg/Point>

#define pressureThreshold   0.5


osg::Geode* drawCube(int cx, int cy, int cz, float size)
{
    osg::Box* unitCube = new osg::Box( osg::Vec3(cx,cy,cz), size);
    osg::ShapeDrawable* unitCubeDrawable = new osg::ShapeDrawable(unitCube);
    osg::Geode* basicShapesGeode = new osg::Geode();
    basicShapesGeode->addDrawable(unitCubeDrawable);

    return basicShapesGeode;
}

osg::Geode* createSphere(float radius, osg::Vec4 color)
{
    osg::Sphere* unitCube = new osg::Sphere(osg::Vec3(0,0,0), radius);
    osg::ShapeDrawable* unitCubeDrawable = new osg::ShapeDrawable(unitCube);

	unitCubeDrawable->setColor(color);

    osg::Geode* basicShapesGeode = new osg::Geode();
    basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON); //Unset shader

    osg::ref_ptr<osg::BlendFunc> blendFunc = new osg::BlendFunc;
    blendFunc->setFunction(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(blendFunc);

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
		_stop = 0;
	}


    // Override the constructor
    virtualHandleCallback(){
        _x = _y = _z = 0;
		virtual_joy_sub_ = nh_.subscribe<std_msgs::Float64MultiArray>("g500/thrusters_input", 1,\
												 &virtualHandleCallback::virtual_joy_Callback, this);
    };
    void operator()(osg::Node* node, osg::NodeVisitor* nv){
        
        osg::Group* currGroup = node->asGroup();
        osg::Node* foundNode;
        
		float size = 0.5;
        
        if (currGroup) {
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
	bool _stop;
	ros::NodeHandle nh_;
	ros::Subscriber virtual_joy_sub_;
};

osg::Group* createVirtualHandle(){
    
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
    upDownCallback(){
	robot_z_sub_ = nh_.subscribe<nav_msgs::Odometry>("uwsim/girona500_odom_RAUVI", 1, &upDownCallback::robotZCallback, this);
    };
	void operator()(osg::Node* node, osg::NodeVisitor* nv){

		osg::Group* currGroup = node->asGroup();
		osg::Node* foundNode;
	
		osg::Geode* geode = dynamic_cast<osg::Geode*>(currGroup->getChild(0));
		osg::Geometry* geometry = geode->getDrawable(0)->asGeometry();

		float a,b;
		if(_state == 0){
			a = 0;
		}
		else{
			a = 1;
		}

		float angle;
		if(_state == -1){
			angle = 90+180;
			b = 0.7;
		}
		else if (_state == 1){
			angle = 90;
			b = 1;
		}

		osg::Vec4Array* colors = new osg::Vec4Array;
		colors->push_back(osg::Vec4(0,0,b,a));

		geometry->setColorArray(colors);
		geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

		osg::MatrixTransform* transform = dynamic_cast<osg::MatrixTransform*>(node);        
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		trans_state = 0;
		_transparency = 1;
		_percentage = 0;
		_state = 0;
		snsrPressure_sub_ = nh_.subscribe<underwater_sensor_msgs::Pressure>("g500/pressure", 1, &pressureWarningCallback::pressureSensorCallback, this);
    };
	void operator()(osg::Node* node, osg::NodeVisitor* nv)
	{
		if(trans_state)
			_transparency += 0.2;
		else
			_transparency -= 0.2;

		if(_transparency >= 1)
			trans_state = 0;
		else if(_transparency <=0)
			trans_state = 1;

		osg::Group* currGroup = node->asGroup();
		osg::Node* foundNode;
	
		osg::Geode* geode = dynamic_cast<osg::Geode*>(currGroup->getChild(0));
		osg::Geometry* geometry = geode->getDrawable(0)->asGeometry();

		float a;
		if(_state == 0)
			a = 0;
		else
			a = _transparency;

		osg::Vec4Array* colors = new osg::Vec4Array;
		colors->push_back(osg::Vec4(1,1,1,a));

		geometry->setColorArray(colors);
		geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

		// Continue callbacks for chlildren too.
		((osg::NodeCallback*) this)->traverse(node, nv);
	};

private:
    float _z, _percentage, _transparency;
    bool trans_state;
    int _state; //-1 -> down // 0 -> deadzone  // 1 -> up
	ros::NodeHandle nh_;
	ros::Subscriber snsrPressure_sub_;
};

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
    text1->setFont("home/jseabra/Downloads/arial.ttf");
    text1->setCharacterSize(characterSize);
    text1->setPosition(osg::Vec3(0,0,0));
    text1->setAxisAlignment(osgText::Text::XZ_PLANE);
    text1->setColor(osg::Vec4(1,0,0,1));
    text1->setText(text);
    text1->setAlignment(osgText::Text::CENTER_CENTER);
//    text1->setDrawMode(osgText::Text::TEXT|osgText::Text::ALIGNMENT|osgText::Text::BOUNDINGBOX);
    text1->setDrawMode(osgText::Text::TEXT);
    //text1->setDataVariance(DYNAMIC); 
	text1->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON);
    geode->addDrawable(text1);

    return geode;
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


class RotationColorCallback : public osg::NodeCallback
{
public:
    // Override the constructor
    RotationColorCallback(){
        _angle = -145;
        state = 0;
    };


    void operator()(osg::Node* node, osg::NodeVisitor* nv){

        if (!state) {
            if (_angle > -145) {
                _angle -= 1;
            }
            else
                state = 1;
        }
        else{
            if (_angle < 150) {
                _angle += 1;
            }
            else
                state = 0;
        }


        osg::MatrixTransform* transform = dynamic_cast<osg::MatrixTransform*>(node);
        
        transform->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(_angle), 0.0, 1.0, 0.0));
        // Continue callbacks for chlildren too.
        ((osg::NodeCallback*) this)->traverse(node, nv);
    };
private:
    bool state;
    float _angle;
    float _r, _g, _b, _a;
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
        _width = width;
        _height = height;
        _percentage = 0;
        state = 0;
    	snsrRange_sub_ = nh_.subscribe<sensor_msgs::Range>("uwsim/g500/range", 1, &sliderGaugeCallback::sensorRngCallback, this);
    };
    
    void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        osg::Group* currGroup = node->asGroup();
        osg::Node* foundNode;
        
        if (currGroup) 
		{
            for (unsigned int i = 0 ; i < currGroup->getNumChildren(); i ++)
            {
                if (currGroup->getChild(i)->getName() == "bar_status")
                {
                    osg::Geode* geode = dynamic_cast<osg::Geode*>(currGroup->getChild(i));
                    osg::Geometry* geometry = geode->getDrawable(0)->asGeometry();
                    
                    float r, g, a;
                    
                    r = _percentage;
                    g = 1.0-_percentage;
                    a = 1.2-_percentage;
                    
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
    bool state;
    float _percentage;
    float _width, _height;

	ros::NodeHandle nh_;
	ros::Subscriber snsrRange_sub_;
};


osg::Geode* createRectangularBox(float width, float height, float line_width, float r, float g, float b, float a){
    
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

osg::Group* createAnalogDial(float up, float right, float size, const std::string& pointer_filename, const std::string& background_filename){
    
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

osg::Group* createBarGauge(float width, float height, float percentage){
    
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



void createWindshield (osg::MatrixTransform *baseTransform)
{
    osg::PositionAttitudeTransform* transform_ = new osg::PositionAttitudeTransform;
    transform_->setPosition(osg::Vec3(0.4, 0, 0.7));
    
    osg::Group* interface = new osg::Group;

    transform_->addChild(interface);

    osg::Group* bar = createBarGauge(0.05, 0.4, 0.2);
    bar->setUpdateCallback(new sliderGaugeCallback(0.05, 0.4));
    transform_->addChild(setNodePosition(0.2, -osg::PI/2-0.4, 0, bar));

    osg::PositionAttitudeTransform* transform_vh = new osg::PositionAttitudeTransform;
    transform_vh->setPosition(osg::Vec3(0.7, 0.5, 0.6));    
    transform_vh->addChild(createVirtualHandle());

    transform_->addChild(transform_vh);

    /*for (float b=0; b>-1*osg::PI/2; b-=10) {
        transform_->addChild(setNodePosition(b, -osg::PI/2+0.5, 0, createTextBox(1 , "windshield")));
    }*/

    osg::MatrixTransform* udArrow =  new osg::MatrixTransform;
    udArrow->addChild(createQuadWithTex(0.2,1,1,1,1,"/home/jcgarsan/.uwsim/data/objects/up.png"));
    udArrow->setUpdateCallback(new upDownCallback());
    osg::MatrixTransform* udArrowTransform = setNodePosition(-0.1, -osg::PI/2+0.4,0,udArrow);
    transform_->addChild(udArrowTransform);

    osg::MatrixTransform* pressureWarning =  new osg::MatrixTransform;
    pressureWarning->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(-90.0), 0.0, 1.0, 0.0));
    pressureWarning->addChild(createQuadWithTex(0.2,1,1,1,1,"/home/jcgarsan/.uwsim/data/objects/warning.png"));
    pressureWarning->setUpdateCallback(new pressureWarningCallback());
    osg::MatrixTransform* pressureWarningTransform = setNodePosition(-0.2, -osg::PI/2,0,pressureWarning);
    transform_->addChild(pressureWarningTransform);

    //This section should be modified to include the user control request icon
/*    osg::MatrixTransform* pressureWarning2 =  new osg::MatrixTransform;
    pressureWarning2->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(-90.0), 0.0, 1.0, 0.0));
    pressureWarning2->addChild(createQuadWithTex(0.2,1,1,1,1,"/home/jcgarsan/.uwsim/data/objects/warning.png"));
    pressureWarning2->setUpdateCallback(new pressureWarningCallback());
    osg::MatrixTransform* pressureWarningTransform2 = setNodePosition(-0.2, -osg::PI/2+0.2, 0, pressureWarning2);
    transform_->addChild(pressureWarningTransform2);*/

    baseTransform->addChild(createBox(osg::Vec3(1.1,0,1.3),0.5,1.5,0.07,osg::Vec4(0.2,0.2,0.2,1)));
    baseTransform->addChild(transform_);
}

#endif
