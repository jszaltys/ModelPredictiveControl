#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

struct Config
{
	Config() : mname()								{}
	Config(const std::string& name) : mname(name)	{}
	void read_cfg(const std::string &filename)
	{
		std::ifstream in;
		in.open(filename.c_str());

		if (in.fail())
		{
			std::ostringstream oss;
			oss << "File: " << filename << " not found\n";
			throw_error(oss.str());
		}

		else
			read_cfg_helper(in);

		in.close();
	}

	std::multimap<std::string, Config>	mcfgmap;
	std::map<std::string, std::string>	mstringmap;
private:
	void read_cfg_helper(std::ifstream &in)
	{
		std::istringstream iss;

		readline(in, iss);
		while (in.good())
		{
			std::string key, name;
			iss >> key;
			if (key == "GROUP")
			{
				if (iss >> name)
				{
					mcfgmap.insert(std::multimap<std::string, Config>::value_type(name, Config(name)))->second.read_cfg_helper(in);
					readline(in, iss);
				}
				else
					throw_error("Missing group name\n");
			}

			else if (key == "END_GROUP")
				break;

			else
			{
				if (key.find_first_not_of(" \t\n\v\f\r") != std::string::npos)
				{
					std::string val;
					iss >> val;
					if (!mstringmap.insert(std::map<std::string, std::string>::value_type(key, val)).second)
					{
						std::ostringstream oss;
						oss << "Key " << key << " already exists in group: " << mname << std::endl;
						throw_error(oss.str());
					}
					readline(in, iss);
				}
				else
					readline(in, iss);
			}
		}
	}
	inline void readline(std::ifstream& in, std::istringstream& iss, const std::string& comm = "//")
	{
		iss.clear();
		std::string buff, tmp;

		std::getline(in, buff, '\n');

		std::size_t pos = buff.find(comm);

		if (pos != std::string::npos)
			tmp = buff.substr(0, pos);
		else
			tmp = buff;

		iss.str(tmp);
	}
	inline void throw_error(const std::string& message)
	{
		throw(message);
	}

	std::string							mname;
};





#endif