#include "master.h"

Master::Master(){}


std::string Master::hasData(const std::string& s)
{
	auto found_null = s.find("null");
	auto b1 		= s.find_first_of("[");
	auto b2 		= s.rfind("}]");
	
	if (found_null != std::string::npos) 
		return "";
	
	else if (b1 != std::string::npos && b2 != std::string::npos) 
		return s.substr(b1, b2 - b1 + 2);
  return "";
}
void Master::run()
{
	h.onMessage			([this](uWS::WebSocket<uWS::SERVER> ws, char *message, size_t length, uWS::OpCode opCode)
	{
		std::string sdata = std::string(message).substr(0, length);
		if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') 
		{	
			auto s = hasData(sdata);	
			if (s != "")
			{
				auto j = nlohmann::json::parse(s);
				
				if (j[0].get<std::string>() == "telemetry")
				{
					nlohmann::json msgJson;

					const std::vector<double> ptsx  = j[1]["ptsx"];
					const std::vector<double> ptsy  = j[1]["ptsy"];
					const double px 				= j[1]["x"];
					const double py 				= j[1]["y"];
					const double psi 				= j[1]["psi"];
					const double v 					= j[1]["speed"]; 
					const double delta 				= j[1]["steering_angle"];
					const double a 					= j[1]["throttle"];
					
					std::vector<double> waypoints_x;
					std::vector<double> waypoints_y;
					std::vector<double> mpc_x_vals;
					std::vector<double> mpc_y_vals;
					std::vector<double> next_x_vals;
					std::vector<double> next_y_vals;
					Eigen::VectorXd	state(mpc.get_state_size());
		
					for (uint32_t i = 0; i < ptsx.size(); i++)
					{
						const double dx = ptsx[i] - px;
						const double dy = ptsy[i] - py;

						waypoints_x.push_back(dx * cos(-psi) - dy * sin(-psi));
						waypoints_y.push_back(dx * sin(-psi) + dy * cos(-psi));
					}
	
					const Eigen::Map<Eigen::VectorXd> x_vals(waypoints_x.data(), waypoints_x.size());
					const Eigen::Map<Eigen::VectorXd> y_vals(waypoints_y.data(), waypoints_y.size());

					Eigen::VectorXd coeffs = polyfit(x_vals, y_vals, mpc.get_polyfitorder());
					
					const double cte			 = polyeval(coeffs, 0);
					const double epsi			 = -std::atan(coeffs[1]);
					
					const double predicted_px 	 = v * mpc.get_dt_constant();
					const double predicted_psi   = v * -delta / mpc.get_Lf_constant() * mpc.get_dt_constant();
					const double predicted_v 	 = v + a * mpc.get_dt_constant();
					const double predicted_cte 	 = cte + v * sin(epsi) * mpc.get_dt_constant();
					const double predicted_epsi  = epsi + v * -delta / mpc.get_Lf_constant() * mpc.get_dt_constant();
					 
					state << predicted_px, 0, predicted_psi, predicted_v, predicted_cte, predicted_epsi;

					
					const std::vector<double> solution = mpc.solve(state, coeffs);
					
					for (uint32_t i = 2; i < solution.size(); i+=2)
					{
						mpc_x_vals.push_back(solution[i]);
						mpc_y_vals.push_back(solution[i + 1]);
					}

					for (uint32_t i = 1; i < 100; ++i)
					{
						next_x_vals.push_back(i);
						next_y_vals.push_back(polyeval(coeffs, i));
					}
			
					msgJson["steering_angle"]	= solution[0] / (deg2rad(25) * mpc.get_Lf_constant());
					msgJson["throttle"]			= solution[1];
					msgJson["mpc_x"]			= mpc_x_vals;
					msgJson["mpc_y"]			= mpc_y_vals;
					msgJson["next_x"]			= next_x_vals;
					msgJson["next_y"]			= next_y_vals;
					
					auto msg = "42[\"steer\"," + msgJson.dump() + "]";
					std::this_thread::sleep_for(std::chrono::milliseconds(100));
					ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
				}
			}
			else
			{
				const std::string msg = "42[\"manual\",{}]";
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			}
		}
	});
	h.onHttpRequest		([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t)
	{
		const std::string s = "<h1>Hello world!</h1>";

		if (req.getUrl().valueLength == 1)
			res->end(s.data(), s.length());

		else
			res->end(nullptr, 0);
	});
	h.onConnection		([this](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req)
	{
		std::cout << "Connected!!!" << std::endl;
	});
	h.onDisconnection	([this](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) 
	{
		ws.close();
		std::cout << "Disconnected" << std::endl;
	});

	if (h.listen(port))
	{
		std::cout << "Listening to port " << port << std::endl;
	}
	else
	{
		std::cerr << "Failed to listen to port" << std::endl;
	}
	h.run();
}

