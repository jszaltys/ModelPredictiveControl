#ifndef __MASTER_H__
#define __MASTER_H__

#include <thread>
#include <chrono>

#include <uWS/uWS.h>
#include "json.hpp"
#include "MPC.h"


class Master
{
public:
	Master();	

	void run					 ();
private:

	std::string hasData			 (const std::string &s);
	uWS::Hub					 h;

	MPC							 mpc;
	
	unsigned int			 	 port = 4567;
};


#endif // __MASTER_H__ 