/*
 * oculuseventhandler.cpp
 *
 *  Created on: Jun 24, 2014
 *      Author: Bjorn Blissing
 */

#include "uwsim/oculuseventhandler.h"

bool OculusEventHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& ad)
{
	switch(ea.getEventType())
	{
	case(osgGA::GUIEventAdapter::KEYDOWN):
		{
			switch(ea.getKey())
			{
			case osgGA::GUIEventAdapter::KEY_R:
				m_oculusViewConfig->resetSensorOrientation();
				return osgGA::GUIEventHandler::handle(ea, ad);
				break;
			default:
				return osgGA::GUIEventHandler::handle(ea, ad);
			} 
		}
	default:
		return osgGA::GUIEventHandler::handle(ea, ad);
	}
}
