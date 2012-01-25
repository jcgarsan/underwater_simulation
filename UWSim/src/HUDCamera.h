/*
 * HUDCamera.h
 * Show a HUD with an image
 *
 *  Created on: 17/05/2011
 *      Author: mprats
 */

#ifndef HUDCAMERA_H_
#define HUDCAMERA_H_

#include "SimulatorConfig.h"
#include "CustomWidget.h"

#include <osg/Image>


/** A ROS Camera subscriber that can be displayed as a widget */
class HUDCamera : public CustomWidget
{
	osgWidget::Widget* widget;

	class widgetUpdateCallback : public osg::Drawable::UpdateCallback {
		osg::Image* image;
		public:
			widgetUpdateCallback(osg::Image *i): osg::Drawable::UpdateCallback() {this->image=i;} 
			virtual void update(osg::NodeVisitor *nv, osg::Drawable *d) {
				osgWidget::Widget * w=static_cast<osgWidget::Widget*>(d);
				w->setImage(image,true,false);			
			}	
	};

public:
	unsigned int width, height;	///< Width and height in pixels of the input image
	unsigned int posx, posy; ///< Default position of the widget, given in pixels wrt bottom-left corner (X to the right, Y upwards)
	double scale;		///< Percentage of default widget scaling (0..1)
	
	osg::Image* osg_image;	//The osg::Image object where to store the ROS image
	bool ready_;	//true if images have been acquired


	/** Constructor from the image and info topics */
	HUDCamera(unsigned int width, unsigned int height, unsigned int posx=0, unsigned int posy=0, double scale=1);

	bool ready() {return ready_;}

	//Interface to be implemented by widgets. Build a widget window with the data to be displayed
	osgWidget::Window* getWidgetWindow();

	~HUDCamera();
};

#endif /* HUDCAMERA_H_ */