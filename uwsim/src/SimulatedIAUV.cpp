/* Copyright (c) 2013 University of Jaume-I.
 * 
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v3.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors:
 *     Mario Prats
 *     Javier Perez
 *     Joao Sebra
 */

#include <uwsim/SceneBuilder.h>
#include <uwsim/osgOceanScene.h>
#include <uwsim/SimulatorConfig.h>
#include <uwsim/SimulatedIAUV.h>
#include <uwsim/URDFRobot.h>
#include <osg/MatrixTransform>
#include <osg/PositionAttitudeTransform>

#include <osg/Point>

#define pressureThreshold	0.5


/** Callback for updating the vehicle lamp according to the vehicle position */
/*
 class LightUpdateCallback:public osg::NodeCallback {
 osg::Transform *trackNode;	///< Node that the light must track

 public:
 LightUpdateCallback(osg::Transform *trackNode)
 {this->trackNode=trackNode;}

 void operator () (osg::Node *node, osg::NodeVisitor *nv) {
 //update light position to track the node
 osg::LightSource *ls=dynamic_cast<osg::LightSource*>(node);
 osg::Light *l=ls->getLight();
 osg::PositionAttitudeTransform *pat=trackNode->asPositionAttitudeTransform();
 osg::Vec3d pos=pat->getPosition();
 l->setPosition( osg::Vec4f(pos.x(),pos.y(),pos.z()-0.5, 1.f) );

 //call to standard callback
 osg::NodeCallback::operator()(node,nv);
 }
 };
 */

/*
 SimulatedIAUV::SimulatedIAUV(osgOcean::OceanScene *oscene, arm_t armtype) {
 vehicle=new SimulatedVehicle(oscene, "GIRONA500/girona500.osg");

 if (armtype==PA10)
 arm=new SimulatedPA10(oscene);
 else if (armtype==ARM5)
 arm=new SimulatedArmFromURDF5(oscene);

 baseTransform=NULL;
 if(vehicle->baseTransform!=NULL && arm->baseTransform!=NULL) {
 baseTransform=vehicle->baseTransform;
 baseTransform->addChild(arm->baseTransform);

 //Vehicle frame to Arm base frame transform
 osg::Matrix m=arm->baseTransform->getMatrix();
 if (armtype==PA10) {
 m.makeRotate(M_PI,1,0,0);
 } else if (armtype==ARM5) {
 }
 arm->baseTransform->setMatrix(m);
 }
 camview=NULL;

 //Set-up a lamp attached to the vehicle
 osg::Light *_light=new osg::Light;
 _light->setLightNum(1);
 _light->setAmbient( osg::Vec4d(1.0f, 1.0f, 1.0f, 1.0f ));
 _light->setDiffuse( osg::Vec4d( 1.0, 1.0, 1.0, 1.0 ) );
 _light->setSpecular(osg::Vec4d( 0.1f, 0.1f, 0.1f, 1.0f ) );
 _light->setDirection(osg::Vec3d(0.0, 0.0, -5.0));
 _light->setSpotCutoff(40.0);
 _light->setSpotExponent(10.0);

 lightSource = new osg::LightSource;
 lightSource->setLight(_light);
 lightSource->setLocalStateSetModes(osg::StateAttribute::ON);
 lightSource->setUpdateCallback(new LightUpdateCallback(baseTransform));
 }
 */

//////////////////////////////////////////////////////////////////////////////

#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/TexMat>
#include <osg/PolygonMode>

#include <osg/NodeVisitor>
#include <osg/Node>

/*
#include <osgDB/ReadFile>

#include <osg/TexGen>
#include <osg/Drawable>
#include <osg/Texture2D>
#include <osg/GLExtensions>
#include <osg/Object>
#include <osg/StateAttribute>
#include <osg/State>
#include <osg/StateSet>
#include <osg/StateAttributeCallback>
#include <osg/Referenced>
#include <osg/BlendFunc>

#include <osg/PositionAttitudeTransform>
#include <osg/Transform>
#include <osg/Matrix>
#include <osg/MatrixTransform>
#include <osg/PositionAttitudeTransform>
#include <osg/TextureRectangle>
#include <osg/Geometry>

#include <osg/NodeCallback>

#include <osgDB/ReadFile>

#include <osgText/Font>
#include <osgText/Text>
*/

osg::Geode* drawCube(int cx, int cy, int cz, float size){

    osg::Box* unitCube = new osg::Box( osg::Vec3(cx,cy,cz), size);
    osg::ShapeDrawable* unitCubeDrawable = new osg::ShapeDrawable(unitCube);
    osg::Geode* basicShapesGeode = new osg::Geode();
    basicShapesGeode->addDrawable(unitCubeDrawable);

    return basicShapesGeode;
}

osg::Geode* createSphere(float radius, osg::Vec4 color){
    osg::Sphere* unitCube = new osg::Sphere(osg::Vec3(0,0,0), radius);
    osg::ShapeDrawable* unitCubeDrawable = new osg::ShapeDrawable(unitCube);

unitCubeDrawable->setColor(color);
//unitCubeDrawable->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::ON);  // Activate blending.

    osg::Geode* basicShapesGeode = new osg::Geode();
    basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON); //Unset shader

    osg::ref_ptr<osg::BlendFunc> blendFunc = new osg::BlendFunc;
    blendFunc->setFunction(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(blendFunc);

    basicShapesGeode->addDrawable(unitCubeDrawable);
    return basicShapesGeode;
}

osg::Geode* createBox(osg::Vec3 center, float tx, float ty, float tz, osg::Vec4 color){

    osg::Box* unitCube = new osg::Box(center, -tx, -ty, -tz);
    osg::ShapeDrawable* unitCubeDrawable = new osg::ShapeDrawable(unitCube);

    unitCubeDrawable->setColor(color);
//unitCubeDrawable->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::ON);  // Activate blending.

    osg::Geode* basicShapesGeode = new osg::Geode();
  //  basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Program(), osg::StateAttribute::ON); //Unset shader

//osg::PolygonMode* polygonMode = new osg::PolygonMode();
//polygonMode->setMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::LINE); 
//basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(polygonMode);

  //  osg::ref_ptr<osg::BlendFunc> blendFunc = new osg::BlendFunc;
   // blendFunc->setFunction(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   // basicShapesGeode->getOrCreateStateSet()->setAttributeAndModes(blendFunc);

    basicShapesGeode->addDrawable(unitCubeDrawable);
    return basicShapesGeode;
}


class virtualHandleCallback : public osg::NodeCallback
{
public:
    void leapCallback(const geometry_msgs::PoseStamped::ConstPtr& posstamped)
	{
		_x = posstamped->pose.position.z;
		_y = posstamped->pose.position.x;
		_z = posstamped->pose.position.y;
		cout << "_x = " << _x << "; _y = " << _y << "; _z = " << _z << endl;

		if(_px == _x && _py == _y && _pz == _z){
			_x = _y = 0; _z = 100;
			_stop = 1;
		}
		else{
			_px = _x;
			_py = _y;
			_pz = _z;
			_stop = 0;
		}		
	} 

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
//		leap_sub_ = nh_.subscribe<geometry_msgs::PoseStamped>("leap_tracker/pose_stamped_out", 1,\
												 &virtualHandleCallback::leapCallback, this);
		virtual_joy_sub_ = nh_.subscribe<std_msgs::Float64MultiArray>("g500/thrusters_input", 1,\
												 &virtualHandleCallback::virtual_joy_Callback, this);
    };
    void operator()(osg::Node* node, osg::NodeVisitor* nv){
        
        osg::Group* currGroup = node->asGroup();
        osg::Node* foundNode;
        
        /*float size = sqrtf(_x*_x + _y*_y + _z*_z)/500;
		if (size < 0.17)
			size = 0.17;*/
		float size = 0.5;
		//cout << "size = " << size << endl;
        
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
			/*
			osg::Group* geode_group = transform->asGroup();
			cout << geode_group->getChild(0)->getName() << endl;
			osg::Geode* geode = dynamic_cast<osg::Geode*>(geode_group->getChild(0));
			osg::ShapeDrawable* drawable = dynamic_cast<osg::ShapeDrawable*>(geode->getDrawable(0));

			if(_stop)
				drawable->setColor(osg::Vec4(0.8,0,0,1));
			else
				drawable->setColor(osg::Vec4(1,1,1,1));	*/		
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
	//ros::Subscriber leap_sub_;
};

osg::Group* createVirtualHandle(){
    
    osg::Cylinder* cylinder = new osg::Cylinder(osg::Vec3(0,0,0.5),0.05,1);
    osg::ShapeDrawable* cylinderDrawable = new osg::ShapeDrawable(cylinder);
    osg::Geode* cylinderGeode = new osg::Geode();
    cylinderGeode->addDrawable(cylinderDrawable);
    cylinderGeode->setName("cylinder");
/*
    osg::Sphere* sphere = new osg::Sphere(osg::Vec3(0,0,0),2);
    osg::ShapeDrawable* sphereDrawable = new osg::ShapeDrawable(sphere);
    osg::Geode* sphereGeode = new osg::Geode();
    sphereGeode->addDrawable(sphereDrawable);
    sphereGeode->setName("sphere");
*/
    osg::PositionAttitudeTransform* transform = new osg::PositionAttitudeTransform;
    transform->addChild(cylinderGeode);
    transform->setName("handle_transform");
/*
	osg::Geode* cube = createBox(osg::Vec3(0,0,0),0.25,0.25,0.1,osg::Vec4(0.5,0.5,0.5,1));
	osg::Geode* cube2 = createBox(osg::Vec3(0,0,0),0.2,0.2,0.15,osg::Vec4(0.7,0.7,0.7,1));
	osg::Geode* cube3 = createBox(osg::Vec3(0,0,0),0.15,0.15,0.2,osg::Vec4(0.9,0.9,0.9,1));
	osg::Geode* cube4 = createBox(osg::Vec3(0,0,0),0.1,0.1,0.25,osg::Vec4(0.95,0.95,0.95,1));
*/
	
	osg::Geode* cube = createSphere(0.15, osg::Vec4(0.8,0.8,0.8,1));

	osg::PositionAttitudeTransform* transform2 = new osg::PositionAttitudeTransform;
	transform2->addChild(cube);
//	transform2->addChild(cube2);
//	transform2->addChild(cube3);
//	transform2->addChild(cube4);
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
		// the movement is inverted!
//		if(abs(_z-odom->pose.pose.position.z) < 0.001)	//Original by Joao
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

		//cout << currGroup->getChild(0)->getName() << endl;
		//currGroup = currGroup->getChild(0)->asGroup();
		//cout << currGroup->getChild(0)->getName() << endl;
		//currGroup = currGroup->getChild(0)->asGroup();
		//cout << currGroup->getChild(0)->getName() << endl;
		//currGroup = currGroup->getChild(0)->asGroup();
		//cout << currGroup->getChild(0)->getName() << endl;
	
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
	    if (abs(pressureValue->pressure) < pressureThreshold)//(pressureValue->pressure > -2)
			_state = 1;
		else
			_state = 0;
		
	}

    // Override the constructor
    pressureWarningCallback(){
	trans_state = 0;
	_transparency = 1;
        _percentage = 0;
        _state = 0;
    	 snsrPressure_sub_ = nh_.subscribe<underwater_sensor_msgs::Pressure>("g500/pressure", 1, &pressureWarningCallback::pressureSensorCallback, this);
    };
	void operator()(osg::Node* node, osg::NodeVisitor* nv){

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

		//cout << currGroup->getChild(0)->getName() << endl;
		//currGroup = currGroup->getChild(0)->asGroup();
		//cout << currGroup->getChild(0)->getName() << endl;
		//currGroup = currGroup->getChild(0)->asGroup();
		//cout << currGroup->getChild(0)->getName() << endl;
		//currGroup = currGroup->getChild(0)->asGroup();
		//cout << currGroup->getChild(0)->getName() << endl;
	
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

osg::Geode* createQuadWithTex(float size, float r, float g, float b, float a, const std::string& filename){

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

osg::Geode* createTextBox(float radius, const std::string& text){

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

osg::MatrixTransform* setNodePosition(float up, float right, float depth, osg::Node* geode){

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
		//cout << rangeValue->range << endl;
	}

    // Override the constructor
    sliderGaugeCallback(float width, float height){
        _width = width;
        _height = height;
        _percentage = 0;
        state = 0;
    	snsrRange_sub_ = nh_.subscribe<sensor_msgs::Range>("uwsim/g500/range", 1, &sliderGaugeCallback::sensorRngCallback, this);
    };
    
    void operator()(osg::Node* node, osg::NodeVisitor* nv){
        
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

/////////////////////////// JOAO
#define SAMPLE_XML_PATH "/home/garciaju/SamplesConfig.xml"

#define CHECK_RC(nRetVal, what)					    \
    if (nRetVal != XN_STATUS_OK)				    \
{								    \
    printf("%s failed: %s\n", what, xnGetStatusString(nRetVal));    \
    return nRetVal;						    \
}

XnBool fileExists(const char *fn)
{
	XnBool exists;
	xnOSDoesFileExist(fn, &exists);
	return exists;
}

Kinect::Kinect(osg::ref_ptr<osg::Vec3Array> vertices, osg::ref_ptr<osg::Vec4Array> colors, unsigned int * pointsI)
{
	float pclRot1 = 1;	//cos(0) = 1
	float pclRot2 = 0;	//sin(0) = 0

	pointsIndices = pointsI;
	pointsXYZ = vertices;
	pointsRGB = colors;
}

Kinect::~Kinect()
{
	g_scriptNode.Release();
	g_DepthGenerator.Release();

	g_ImageGenerator.Release();
	g_Context.Release();
}

int Kinect::init()
{
	XnStatus nRetVal = XN_STATUS_OK;
	xn::EnumerationErrors errors;

	const char *fn = NULL;
	if    (fileExists(SAMPLE_XML_PATH)) fn = SAMPLE_XML_PATH;
	else {
	printf("Could not find '%s'. Aborting.\n" , SAMPLE_XML_PATH);
	return XN_STATUS_ERROR;
	}
	printf("Reading config from: '%s'\n", fn);

	nRetVal = g_Context.InitFromXmlFile(fn, g_scriptNode, &errors);
	if (nRetVal == XN_STATUS_NO_NODE_PRESENT)
	{
		XnChar strError[1024];
		errors.ToString(strError, 1024);
		printf("%s\n", strError);
		return (nRetVal);
	}
	else if (nRetVal != XN_STATUS_OK)
	{
		printf("Open failed: %s\n", xnGetStatusString(nRetVal));
		return (nRetVal);
	}

	nRetVal = g_Context.FindExistingNode(XN_NODE_TYPE_DEPTH, g_DepthGenerator);
	CHECK_RC(nRetVal,"No depth");

	nRetVal = g_Context.FindExistingNode(XN_NODE_TYPE_IMAGE, g_ImageGenerator);
	CHECK_RC(nRetVal, "No color");

	//calibrate depth camera to match RGB camera
	xnSetViewPoint(g_DepthGenerator,g_ImageGenerator);

	nRetVal = g_Context.StartGeneratingAll();
	CHECK_RC(nRetVal, "StartGenerating");
}

void Kinect::update()
{
	//get the updated data from kinect
	depthData = g_DepthGenerator.GetDepthMap();
	imageData = g_ImageGenerator.GetRGB24ImageMap();
	
	pointsXYZ->clear();
	pointsRGB->clear();

	int x = 0;
	int z = 0;

	//this is to rotate the entire point cloud
	pclRot1 = cos(DEG2RAD*kinectAngle);
	pclRot2 = sin(DEG2RAD*kinectAngle);

	//cycle through all points, insert them into the RGB array and insert them into the projective coordinates array (points2D)
	for(int i=0; i<640*480; i++)
	{
		x++;
		if(x == 640)
		{
			z++;
			x=0;
		}

		//get depth data, place it on image plane and convert to real world coordinates
		points2D[i].X = x;
		points2D[i].Y = z;
		points2D[i].Z = depthData[i];

		//get RGB data
		pointsRGB->push_back(osg::Vec4f(imageData[i].nRed/255.0, imageData[i].nGreen/255.0, imageData[i].nBlue/255.0,1.0f));

		pointsIndices[i] = i; //add this point to list of points to draw
	}
 	
	//convert the whole projective coordinates array into real world coordinates (points3D)
	g_DepthGenerator.ConvertProjectiveToRealWorld(307200, points2D, points3D);

	//apply kinect->openAR transform and rotate all points matching kinect/table angle
	for(int i=0; i<640*480; i++)
	{
		pointsXYZ->push_back( osg::Vec3f(points3D[i].X/1000.0, (-points3D[i].Z*pclRot1 - points3D[i].Y*pclRot2)/1000.0, (-points3D[i].Z*pclRot2 + points3D[i].Y*pclRot1)/1000.0));
	}
}

//this function performs a 3DHough transform to detect the angle that the kinect is making with the scene's
//biggest plane (the table). It then returns and object with the plane's spherical coordinates.
//The value of Phi is then used to define "kinectAngle" which is used to rotate the point cloud
/*TableTransform Kinect::tableDetection()
{
	float tolerance = 0.001;
	
	float thetaStep = 0.1;
	float phiStep = 0.005;
	float rohStep = 2;
	float rohRange = 600;
	float rohStart = 400;

	int nTheta = (2*M_PI)/thetaStep;
	int nPhi = M_PI/phiStep;
	int nRoh = rohRange/rohStep;

	printf("nTheta:%d nPhi:%d nRoh:%d\n", nTheta, nPhi, nRoh);

	float theta[nTheta];
	float phi[nPhi];
	float roh[nRoh];

	//  Allocate 3D Array
	int ***accGrid = new int**[nTheta];
	for(int i = 0; i < nTheta; i++)
	{
		accGrid[i] = new int*[nPhi];

		for(int j = 0; j < nPhi; j++)
		{
		    accGrid[i][j] = new int[nRoh];
		}
	}

	for(int i=0; i<nTheta; i++)
	{
		theta[i] = i*thetaStep;
	}

	for(int i=0; i<nPhi; i++)
	{
		phi[i] = i*phiStep;
	}

	for(int i=0; i<nRoh; i++)
	{
		roh[i] = i*rohStep+rohStart;
	}

	for(int j=0; j<nTheta; j++)
	{
		for(int k=0; k<nPhi; k++)
		{
			for(int l=0; l<nRoh; l++)
			{
				accGrid[j][k][l] = 0;					
			}
		}
	}

	float calculatedRoh;
	float result;

	for(int i=0; i<640*480; i+=400)
	{
		for(int j=0; j<nTheta; j++)
		{
			for(int k=0; k<nPhi; k++)
			{
				calculatedRoh = pointsXYZ[i].x*cos(theta[j])*sin(phi[k])+pointsXYZ[i].y*sin(phi[k])*sin(theta[j])+pointsXYZ[i].z*cos(phi[k]);

				for(int l=0; l<nRoh; l++)
				{
					result = (calculatedRoh)/roh[l];

					if(result <= 1+tolerance && result >= 1-tolerance)
					{
						accGrid[j][k][l] ++;						
					}
				}
			}
		}
	}

	int max = 0;
	int maxJ,maxK,maxL;

	for(int j=0; j<nTheta; j++)
	{
		for(int k=0; k<nPhi; k++)
		{
			for(int l=0; l<nRoh; l++)
			{
				if(accGrid[j][k][l] > max)
				{
					max = accGrid[j][k][l];
					maxJ = j;
					maxK = k;
					maxL = l;
				}
			}
		}
	}

	printf("\n\ntheta: %f phi: %f roh: %f \nj: %d k: %d l: %d max: %d \n\n",theta[maxJ], phi[maxK], roh[maxL], maxJ, maxK, maxL, max);
	printf("theta(degrees): %f phi(degrees): %f\n",theta[maxJ]/DEG2RAD, phi[maxJ]/DEG2RAD);

	//deallocate 3d array
	for(int i = 0; i < nTheta; i++)
	{
		for(int j = 0; j < nPhi; j++)
		{
			delete[] accGrid[i][j];
		}

		delete[] accGrid[i];
	}
	delete[] accGrid;

	table.phi = phi[maxK]/DEG2RAD;
	table.theta = theta[maxJ]/DEG2RAD;
	table.roh = roh[maxL];

	kinectAngle = 180-table.phi;

	return table;
}*/
/////////////////////////

//////////////////////////////////////////////////////////////////////////////

SimulatedIAUV::SimulatedIAUV(SceneBuilder *oscene, Vehicle vehicleChars) : urdf(new URDFRobot(oscene->scene->getOceanScene(), vehicleChars))
{
  name = vehicleChars.name;
  baseTransform = new osg::MatrixTransform;

  if (urdf->baseTransform != NULL /* && arm->baseTransform!=NULL*/)
  {
    baseTransform->addChild(urdf->baseTransform);
    baseTransform->setName(vehicleChars.name);
  }
  scale=osg::Vec3d(vehicleChars.scale[0],vehicleChars.scale[1],vehicleChars.scale[2]);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (1)
    {
	osg::PositionAttitudeTransform* transform_ = new osg::PositionAttitudeTransform;
	transform_->setPosition(osg::Vec3(0.4, 0, 0.7));
	
	osg::Group* interface = new osg::Group;

	transform_->addChild(interface);

	osg::Group* bar = createBarGauge(0.05, 0.4, 0.2);
	bar->setUpdateCallback(new sliderGaugeCallback(0.05, 0.4));
	transform_->addChild(setNodePosition(0.2, -osg::PI/2-0.4, 0, bar));

	//transform_->addChild(createAnalogDial(0.2, -osg::PI/2, 0.1, "/home/jseabra/Downloads/2.png", "/home/jseabra/Downloads/bottom.png"));

	osg::PositionAttitudeTransform* transform_vh = new osg::PositionAttitudeTransform;
	transform_vh->setPosition(osg::Vec3(0.7, 0.5, 0.6));	
	transform_vh->addChild(createVirtualHandle());

	transform_->addChild(transform_vh);

	/*for (float b=0; b>-1*osg::PI/2; b-=10) {
		transform_->addChild(setNodePosition(b, -osg::PI/2+0.5, 0, createTextBox(1 , "windshield")));
	}*/

	osg::MatrixTransform* udArrow =  new osg::MatrixTransform;
//	udArrow->addChild(createQuadWithTex(0.2,1,1,1,1,"/home/jseabra/Downloads/up.png"));
	udArrow->addChild(createQuadWithTex(0.2,1,1,1,1,"/home/jcgarsan/.uwsim/data/objects/up.png"));
	udArrow->setUpdateCallback(new upDownCallback());
	osg::MatrixTransform* udArrowTransform = setNodePosition(-0.1, -osg::PI/2+0.4,0,udArrow);
	transform_->addChild(udArrowTransform);

	osg::MatrixTransform* pressureWarning =  new osg::MatrixTransform;
	pressureWarning->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(-90.0), 0.0, 1.0, 0.0));
//	pressureWarning->addChild(createQuadWithTex(0.2,1,1,1,1,"/home/jseabra/Downloads/warning.png"));
	pressureWarning->addChild(createQuadWithTex(0.2,1,1,1,1,"/home/jcgarsan/.uwsim/data/objects/warning.png"));
	pressureWarning->setUpdateCallback(new pressureWarningCallback());
	osg::MatrixTransform* pressureWarningTransform = setNodePosition(-0.2, -osg::PI/2,0,pressureWarning);
	transform_->addChild(pressureWarningTransform);


	baseTransform->addChild(createBox(osg::Vec3(1.1,0,1.3),0.5,1.5,0.07,osg::Vec4(0.2,0.2,0.2,1)));
//	baseTransform->addChild(createBox(osg::Vec3(0,0,0),3,3,3,osg::Vec4(1,1,1,0.1)));
	baseTransform->addChild(transform_);
	
	vertices = osg::ref_ptr < osg::Vec3Array > (new osg::Vec3Array());
	colors = osg::ref_ptr < osg::Vec4Array > (new osg::Vec4Array());
	
/*    kinect= new Kinect(vertices,colors,pointsIndices);
	kinect->init();
	kinect->update();
	
	geode = osg::ref_ptr < osg::Geode > (new osg::Geode());
    geometry = osg::ref_ptr < osg::Geometry > (new osg::Geometry());

    geometry->setVertexArray(vertices);
    geometry->setColorArray(colors);
    geometry->setColorBinding(osg::Geometry::BIND_PER_VERTEX);

    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, 307200));

    geode->addDrawable(geometry.get());
    osg::StateSet* state = geometry->getOrCreateStateSet();
    state->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

    osg::Point *point = new osg::Point;
    point->setSize(2);
    state->setAttribute(point);
           
    geode->setUpdateCallback(new KinectCallback(kinect,geometry));
    
   	osg::MatrixTransform* kinectTransform =  new osg::MatrixTransform();
   	osg::Matrixd  kinectMatrix ;
   	
   	kinectMatrix.makeRotate(osg::Quat(3.14, osg::Vec3d(1, 0, 0), 0, osg::Vec3d(0, 1, 0), -1.57,osg::Vec3d(0, 0, 1)));

    kinectTransform->setMatrix(kinectMatrix);
    kinectTransform->addChild(geode);
    
    geode->setNodeMask(0x04 | 0x02 | 0x01);
    
    baseTransform->addChild(kinectTransform);*/

	

	

	
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //Add virtual  cameras in config file
  while (vehicleChars.Vcams.size() > 0)
  {
    Vcam vcam = vehicleChars.Vcams.front();
    OSG_INFO << "Adding a virtual camera " << vcam.name << "..." << std::endl;
    vehicleChars.Vcams.pop_front();
    //Camera frame given wrt vehicle origin frame.
    //Remember that in opengl/osg, the camera frame is a right-handed system with Z going backwards (opposite to the viewing direction) and Y up.
    osg::ref_ptr < osg::Transform > vMc = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMc->asPositionAttitudeTransform()->setPosition(osg::Vec3d(vcam.position[0], vcam.position[1], vcam.position[2]));
    vMc->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(vcam.orientation[0], osg::Vec3d(1, 0, 0), vcam.orientation[1], osg::Vec3d(0, 1, 0),
                  vcam.orientation[2], osg::Vec3d(0, 0, 1)));
    urdf->link[vcam.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMc);
    camview.push_back(
        VirtualCamera(oscene->root, vcam.name, vcam.linkName, vMc, vcam.resw, vcam.resh, vcam.baseLine, vcam.frameId,
                      vcam.fov,oscene,vcam.std,vcam.parameters.get(), 0, vcam.bw));
    if (vcam.showpath)
      camview[camview.size() - 1].showPath(vcam.showpath);
    OSG_INFO << "Done adding a virtual camera..." << std::endl;
  }

  //Add virtual range cameras in config file
  while (vehicleChars.VRangecams.size() > 0)
  {
    Vcam vcam = vehicleChars.VRangecams.front();
    OSG_INFO << "Adding a virtual camera " << vcam.name << "..." << std::endl;
    vehicleChars.VRangecams.pop_front();
    //Camera frame given wrt vehicle origin frame.
    //Remember that in opengl/osg, the camera frame is a right-handed system with Z going backwards (opposite to the viewing direction) and Y up.
    osg::ref_ptr < osg::Transform > vMc = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMc->asPositionAttitudeTransform()->setPosition(osg::Vec3d(vcam.position[0], vcam.position[1], vcam.position[2]));
    vMc->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(vcam.orientation[0], osg::Vec3d(1, 0, 0), vcam.orientation[1], osg::Vec3d(0, 1, 0),
                  vcam.orientation[2], osg::Vec3d(0, 0, 1)));
    urdf->link[vcam.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMc);
    camview.push_back(
        VirtualCamera(oscene->root, vcam.name, vcam.linkName, vMc, vcam.resw, vcam.resh, vcam.baseLine, vcam.frameId,
                      vcam.fov,NULL,0,vcam.parameters.get(), 1, 0));
    if (vcam.showpath)
      camview[camview.size() - 1].showPath(vcam.showpath);
    //Check underwaterParticles to change the mask
    if (!vcam.underwaterParticles)
      camview[camview.size() - 1].textureCamera->setCullMask(oscene->scene->getOceanScene()->getNormalSceneMask());
    OSG_INFO << "Done adding a virtual camera..." << std::endl;
  }

  // Adding Structured light projector
  while (vehicleChars.sls_projectors.size() > 0)
  {
    OSG_INFO << "Adding a structured light projector..." << std::endl;
    slProjector slp;
    slp = vehicleChars.sls_projectors.front();
    vehicleChars.sls_projectors.pop_front();
    osg::ref_ptr < osg::Transform > vMp = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMp->asPositionAttitudeTransform()->setPosition(osg::Vec3d(slp.position[0], slp.position[1], slp.position[2]));
    vMp->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(slp.orientation[0], osg::Vec3d(1, 0, 0), slp.orientation[1], osg::Vec3d(0, 1, 0), slp.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    urdf->link[slp.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMp);
    //camview.push_back(VirtualCamera(oscene->root, "slp_camera", vMp, 512, 512,slp.fov,102.4));
    sls_projectors.push_back(VirtualSLSProjector(slp.name, slp.linkName, oscene->root, //maybe oscene->scene->localizedWorld ?
                                                 vMp, slp.image_name, slp.fov, (slp.laser) ? true : false));
    camview.push_back(sls_projectors.back().camera);
    OSG_INFO << "Done adding a structured light projector..." << std::endl;
  }

  //Adding range sensors
  while (vehicleChars.range_sensors.size() > 0)
  {
    OSG_INFO << "Adding a virtual range sensor..." << std::endl;
    rangeSensor rs;
    rs = vehicleChars.range_sensors.front();
    vehicleChars.range_sensors.pop_front();
    osg::ref_ptr < osg::Transform > vMr = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMr->asPositionAttitudeTransform()->setPosition(osg::Vec3d(rs.position[0], rs.position[1], rs.position[2]));
    vMr->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(rs.orientation[0], osg::Vec3d(1, 0, 0), rs.orientation[1], osg::Vec3d(0, 1, 0), rs.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    urdf->link[rs.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMr);
    range_sensors.push_back(
        VirtualRangeSensor(rs.name, rs.linkName, oscene->scene->localizedWorld, vMr, rs.range, (rs.visible) ? true : false,
        oscene->scene->getOceanScene()->getARMask()));
    OSG_INFO << "Done adding a virtual range sensor..." << std::endl;
  }

  //Adding imus
  while (vehicleChars.imus.size() > 0)
  {
    OSG_INFO << "Adding an IMU..." << std::endl;
    Imu imu;
    imu = vehicleChars.imus.front();
    vehicleChars.imus.pop_front();
    osg::ref_ptr < osg::Transform > vMi = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMi->asPositionAttitudeTransform()->setPosition(osg::Vec3d(imu.position[0], imu.position[1], imu.position[2]));
    vMi->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(imu.orientation[0], osg::Vec3d(1, 0, 0), imu.orientation[1], osg::Vec3d(0, 1, 0), imu.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    urdf->link[imu.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMi);
    imus.push_back(InertialMeasurementUnit(imu.name, imu.linkName, vMi, oscene->scene->localizedWorld->getMatrix(), imu.std));
    OSG_INFO << "Done adding an IMU..." << std::endl;
  }

  //Adding pressure sensors
  while (vehicleChars.pressure_sensors.size() > 0)
  {
    OSG_INFO << "Adding a pressure sensor..." << std::endl;
    XMLPressureSensor ps;
    ps = vehicleChars.pressure_sensors.front();
    vehicleChars.pressure_sensors.pop_front();
    osg::ref_ptr < osg::Transform > vMs = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMs->asPositionAttitudeTransform()->setPosition(osg::Vec3d(ps.position[0], ps.position[1], ps.position[2]));
    vMs->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(ps.orientation[0], osg::Vec3d(1, 0, 0), ps.orientation[1], osg::Vec3d(0, 1, 0), ps.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    urdf->link[ps.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMs);
    pressure_sensors.push_back(PressureSensor(ps.name, ps.linkName, vMs, oscene->scene->localizedWorld->getMatrix(), ps.std));
    OSG_INFO << "Done adding an Pressure Sensor..." << std::endl;
  }

  //Adding GPS sensors
  while (vehicleChars.gps_sensors.size() > 0)
  {
    OSG_INFO << "Adding a gps sensor..." << std::endl;
    XMLGPSSensor ps;
    ps = vehicleChars.gps_sensors.front();
    vehicleChars.gps_sensors.pop_front();
    osg::ref_ptr < osg::Transform > vMs = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMs->asPositionAttitudeTransform()->setPosition(osg::Vec3d(ps.position[0], ps.position[1], ps.position[2]));
    vMs->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(ps.orientation[0], osg::Vec3d(1, 0, 0), ps.orientation[1], osg::Vec3d(0, 1, 0), ps.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    urdf->link[ps.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMs);
    gps_sensors.push_back(GPSSensor(oscene->scene, ps.name, ps.linkName , vMs, oscene->scene->localizedWorld->getMatrix(), ps.std));
    OSG_INFO << "Done adding an GPS Sensor..." << std::endl;
  }

  //Adding dvl sensors
  while (vehicleChars.dvl_sensors.size() > 0)
  {
    OSG_INFO << "Adding a dvl sensor..." << std::endl;
    XMLDVLSensor ps;
    ps = vehicleChars.dvl_sensors.front();
    vehicleChars.dvl_sensors.pop_front();
    osg::ref_ptr < osg::Transform > vMs = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMs->asPositionAttitudeTransform()->setPosition(osg::Vec3d(ps.position[0], ps.position[1], ps.position[2]));
    vMs->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(ps.orientation[0], osg::Vec3d(1, 0, 0), ps.orientation[1], osg::Vec3d(0, 1, 0), ps.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    urdf->link[ps.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMs);
    dvl_sensors.push_back(DVLSensor(ps.name, ps.linkName, vMs, oscene->scene->localizedWorld->getMatrix(), ps.std));
    OSG_INFO << "Done adding an DVL Sensor..." << std::endl;
  }

  //Adding Multibeam sensors
  while (vehicleChars.multibeam_sensors.size() > 0)
  {
    OSG_INFO << "Adding a Multibeam sensor..." << std::endl;
    XMLMultibeamSensor MB;
    MB = vehicleChars.multibeam_sensors.front();
    vehicleChars.multibeam_sensors.pop_front();
    osg::ref_ptr < osg::Transform > vMs = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMs->asPositionAttitudeTransform()->setPosition(osg::Vec3d(MB.position[0], MB.position[1], MB.position[2]));
    vMs->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(MB.orientation[0], osg::Vec3d(1, 0, 0), MB.orientation[1], osg::Vec3d(0, 1, 0), MB.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    urdf->link[MB.link]->getParent(0)->getParent(0)->asGroup()->addChild(vMs);
    unsigned int mask;
    if(MB.underwaterParticles)
      mask=oscene->scene->getOceanScene()->getARMask();
    else
      mask=oscene->scene->getOceanScene()->getNormalSceneMask(); //Normal Scene mask should be enough for range sensor
    MultibeamSensor mb = MultibeamSensor(oscene->root, MB.name, MB.linkName, vMs, MB.initAngle, MB.finalAngle, MB.angleIncr,
                                         MB.range,mask,MB.visible,mask);
    multibeam_sensors.push_back(mb);
    for(unsigned int i=0;i<mb.nCams;i++)
      camview.push_back(mb.vcams[i]);
    OSG_INFO << "Done adding a Multibeam Sensor..." << std::endl;
  }

  //Adding object pickers
  while (vehicleChars.object_pickers.size() > 0)
  {
    OSG_INFO << "Adding an object picker..." << std::endl;
    rangeSensor rs;
    rs = vehicleChars.object_pickers.front();
    vehicleChars.object_pickers.pop_front();
    osg::ref_ptr < osg::Transform > vMr = (osg::Transform*)new osg::PositionAttitudeTransform;
    vMr->asPositionAttitudeTransform()->setPosition(osg::Vec3d(rs.position[0], rs.position[1], rs.position[2]));
    vMr->asPositionAttitudeTransform()->setAttitude(
        osg::Quat(rs.orientation[0], osg::Vec3d(1, 0, 0), rs.orientation[1], osg::Vec3d(0, 1, 0), rs.orientation[2],
                  osg::Vec3d(0, 0, 1)));
    vMr->setName("ObjectPickerNode");
    urdf->link[rs.link]->asGroup()->addChild(vMr);
    object_pickers.push_back(ObjectPicker(rs.name, oscene->scene->localizedWorld, vMr, rs.range, true, urdf,
        oscene->scene->getOceanScene()->getARMask()));
    OSG_INFO << "Done adding an object picker..." << std::endl;
  }

  devices.reset(new SimulatedDevices());
  devices->applyConfig(this, vehicleChars, oscene);

  //Set-up a lamp attached to the vehicle: TODO
  /*
   osg::Light *_light=new osg::Light;
   _light->setLightNum(1);
   _light->setAmbient( osg::Vec4d(1.0f, 1.0f, 1.0f, 1.0f ));
   _light->setDiffuse( osg::Vec4d( 1.0, 1.0, 1.0, 1.0 ) );
   _light->setSpecular(osg::Vec4d( 0.1f, 0.1f, 0.1f, 1.0f ) );
   _light->setDirection(osg::Vec3d(0.0, 0.0, -5.0));
   _light->setSpotCutoff(40.0);
   _light->setSpotExponent(10.0);

   lightSource = new osg::LightSource;
   lightSource->setLight(_light);
   lightSource->setLocalStateSetModes(osg::StateAttribute::ON);
   lightSource->setUpdateCallback(new LightUpdateCallback(baseTransform));
   */
}

/*
 void SimulatedIAUV::setVirtualCamera(std::string name, osg::Transform* transform, int width, int height) {
 //Set I-AUV virtual camera

 baseTransform->asGroup()->addChild(transform);

 if (camview==NULL)
 camview=new VirtualCamera(name, transform, width, height);
 }
 */

/** Sets the vehicle position. (x,y,z) given wrt to the world frame. (roll,pitch,yaw) are RPY angles in the local frame */
void SimulatedIAUV::setVehiclePosition(double x, double y, double z, double roll, double pitch, double yaw)
{
  osg::Matrixd S, T, Rx, Ry, Rz, transform;
  T.makeTranslate(x, y, z);
  Rx.makeRotate(roll, 1, 0, 0);
  Ry.makeRotate(pitch, 0, 1, 0);
  Rz.makeRotate(yaw, 0, 0, 1);
  S.makeScale(scale);
  transform = S * Rz * Ry * Rx * T;
  setVehiclePosition(transform);
}

void SimulatedIAUV::setVehiclePosition(osg::Matrixd m)
{
  baseTransform->setMatrix(m);
}

