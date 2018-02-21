#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "PID.h"

using namespace std;

using namespace Eigen;
typedef Matrix<double, Dynamic, 7, RowMajor > RowMatrixXi;


// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  } 
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}



int main() {
  uWS::Hub h;




  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // Start in Lane # 1
  int lane = 1;        

  // velocity set point       
  double ref_vel = .0; // mph

  h.onMessage([&lane,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&ref_vel](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

                double car_b_speed = .0;
                double car_b_s = 0.0;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	
                vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];


/*              
//                RowMatrixXi sefu; 




                MatrixXd sefu(sensor_fusion.size(),7);


                sefu << sensor_fusion[0][0], sensor_fusion[0][1], sensor_fusion[0][2], sensor_fusion[0][3], sensor_fusion[0][4], sensor_fusion[0][5], sensor_fusion[0][6],
                        sensor_fusion[1][0], sensor_fusion[1][1], sensor_fusion[1][2], sensor_fusion[1][3], sensor_fusion[1][4], sensor_fusion[1][5], sensor_fusion[1][6],
                        sensor_fusion[2][0], sensor_fusion[2][1], sensor_fusion[2][2], sensor_fusion[2][3], sensor_fusion[2][4], sensor_fusion[2][5], sensor_fusion[2][6],
                        sensor_fusion[3][0], sensor_fusion[3][1], sensor_fusion[3][2], sensor_fusion[3][3], sensor_fusion[3][4], sensor_fusion[3][5], sensor_fusion[3][6],
                        sensor_fusion[4][0], sensor_fusion[4][1], sensor_fusion[4][2], sensor_fusion[4][3], sensor_fusion[4][4], sensor_fusion[4][5], sensor_fusion[4][6],
                        sensor_fusion[5][0], sensor_fusion[5][1], sensor_fusion[5][2], sensor_fusion[5][3], sensor_fusion[5][4], sensor_fusion[5][5], sensor_fusion[5][6],
                        sensor_fusion[6][0], sensor_fusion[6][1], sensor_fusion[6][2], sensor_fusion[6][3], sensor_fusion[6][4], sensor_fusion[6][5], sensor_fusion[6][6],
                        sensor_fusion[7][0], sensor_fusion[7][1], sensor_fusion[7][2], sensor_fusion[7][3], sensor_fusion[7][4], sensor_fusion[7][5], sensor_fusion[7][6],
                        sensor_fusion[8][0], sensor_fusion[8][1], sensor_fusion[8][2], sensor_fusion[8][3], sensor_fusion[8][4], sensor_fusion[8][5], sensor_fusion[8][6],
                        sensor_fusion[9][0], sensor_fusion[9][1], sensor_fusion[9][2], sensor_fusion[9][3], sensor_fusion[9][4], sensor_fusion[9][5], sensor_fusion[9][6],
                        sensor_fusion[10][0], sensor_fusion[10][1], sensor_fusion[10][2], sensor_fusion[10][3], sensor_fusion[10][4], sensor_fusion[10][5], sensor_fusion[10][6],
                        sensor_fusion[11][0], sensor_fusion[11][1], sensor_fusion[11][2], sensor_fusion[11][3], sensor_fusion[11][4], sensor_fusion[11][5], sensor_fusion[11][6];
 


                std::ptrdiff_t i, j;

                double thisone= sefu.minCoeff(&i, &j);
                  
                VectorXd this_vector(7);

                this_vector = sefu.row(i);  
               
                std::cout << this_vector << std::endl;

*/              

                //auto Veh_ahead(sensor_fusion.size());
                //std::vector<vector<double>> Veh_ahead(sensor_fusion.size());
                

                // Get sefu_VehAhead from sensor_fusion
                // a vector of vectors with all vehicles at 300 meters ahead

                vector<vector<double>> sefu_VehAhead = j[1]["sensor_fusion"];
                     
                int p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {
                    
                    sefu_VehAhead[i][0] = 0.0;
                    sefu_VehAhead[i][1] = 0.0;
                    sefu_VehAhead[i][2] = 0.0;
                    sefu_VehAhead[i][3] = 0.0;
                    sefu_VehAhead[i][4] = 0.0;
                    sefu_VehAhead[i][5] = 0.0;
                    sefu_VehAhead[i][6] = 0.0;
                    
                    if (((sensor_fusion[i][5] - car_s) < 400.0) && ((sensor_fusion[i][5] - car_s) > 0))
                        {
                                         
                        sefu_VehAhead[p][0] = sensor_fusion[i][0];
                        sefu_VehAhead[p][1] = sensor_fusion[i][1];
                        sefu_VehAhead[p][2] = sensor_fusion[i][2];
                        sefu_VehAhead[p][3] = sensor_fusion[i][3];
                        sefu_VehAhead[p][4] = sensor_fusion[i][4];
                        sefu_VehAhead[p][5] = sensor_fusion[i][5];
                        sefu_VehAhead[p][6] = sensor_fusion[i][6];
                        p++;
                        
                        }
                }

                
                // Get VehAhead_same from sefu_VehAhead
                // a vector of vectors with all vehicles at 300 meters ahead and on same lane

                vector<vector<double>> VehAhead_SameLane = j[1]["sensor_fusion"];
                     
                p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {
                    
                    VehAhead_SameLane[i][0] = 0.0;
                    VehAhead_SameLane[i][1] = 0.0;
                    VehAhead_SameLane[i][2] = 0.0;
                    VehAhead_SameLane[i][3] = 0.0;
                    VehAhead_SameLane[i][4] = 0.0;
                    VehAhead_SameLane[i][5] = 0.0;
                    VehAhead_SameLane[i][6] = 0.0;
                    
                    if ((abs(sefu_VehAhead[i][6] - car_d) < 2.0))
                        {
                                         
                        VehAhead_SameLane[p][0] = sefu_VehAhead[i][0];
                        VehAhead_SameLane[p][1] = sefu_VehAhead[i][1];
                        VehAhead_SameLane[p][2] = sefu_VehAhead[i][2];
                        VehAhead_SameLane[p][3] = sefu_VehAhead[i][3];
                        VehAhead_SameLane[p][4] = sefu_VehAhead[i][4];
                        VehAhead_SameLane[p][5] = sefu_VehAhead[i][5];
                        VehAhead_SameLane[p][6] = sefu_VehAhead[i][6];
                        p++;
                         
                        }
                }

               
 
                // Get SameLane_abc from VehAhead_SameLane
                // a vector with telemetry of the first 3 vehicle in order of appereance 
 

                vector<vector<double>> SameLane_abc(3,vector<double>(7));

 
                if((VehAhead_SameLane[1][5] == 0) && (VehAhead_SameLane[2][5] == 0))
                    for(int i=0; i < 7 ; i++) SameLane_abc[0][i] = VehAhead_SameLane[0][i];   
                else  if(VehAhead_SameLane[2][5] == 0)
                           {    
                           if(VehAhead_SameLane[0][5] < VehAhead_SameLane[1][5])
                               {
                               for(int i=0; i < 7 ; i++) SameLane_abc[0][i] = VehAhead_SameLane[0][i]; 
                               for(int i=0; i < 7 ; i++) SameLane_abc[1][i] = VehAhead_SameLane[1][i];
                               }   
                           else
                               {  
                               for(int i=0; i < 7 ; i++) SameLane_abc[0][i] = VehAhead_SameLane[1][i];
                               for(int i=0; i < 7 ; i++) SameLane_abc[1][i] = VehAhead_SameLane[0][i]; 
                               }  
                           }    
                       
                else { 


                      if((VehAhead_SameLane[0][5] < VehAhead_SameLane[1][5])  && 
                         (VehAhead_SameLane[0][5] < VehAhead_SameLane[2][5]))   
                         for(int i=0; i < 7 ; i++) SameLane_abc[0][i] = VehAhead_SameLane[0][i];   

                      else if((VehAhead_SameLane[0][5] > VehAhead_SameLane[1][5]) &&
                              (VehAhead_SameLane[0][5] > VehAhead_SameLane[2][5]))  
                              for(int i=0; i < 7 ; i++) SameLane_abc[2][i] = VehAhead_SameLane[0][i];   

                      else for(int i=0; i < 7 ; i++) SameLane_abc[1][i] = VehAhead_SameLane[0][i];   
               

                      if((VehAhead_SameLane[1][5] < VehAhead_SameLane[0][5])  && 
                         (VehAhead_SameLane[1][5] < VehAhead_SameLane[2][5]))   
                         for(int i=0; i < 7 ; i++) SameLane_abc[0][i] = VehAhead_SameLane[1][i];   

                      else if((VehAhead_SameLane[1][5] > VehAhead_SameLane[0][5]) &&
                              (VehAhead_SameLane[1][5] > VehAhead_SameLane[2][5]))  
                              for(int i=0; i < 7 ; i++) SameLane_abc[2][i] = VehAhead_SameLane[1][i];   

                      else for(int i=0; i < 7 ; i++) SameLane_abc[1][i] = VehAhead_SameLane[1][i];   


                      if((VehAhead_SameLane[2][5] < VehAhead_SameLane[0][5])  && 
                         (VehAhead_SameLane[2][5] < VehAhead_SameLane[1][5]))   
                         for(int i=0; i < 7 ; i++) SameLane_abc[0][i] = VehAhead_SameLane[2][i];   

                      else if((VehAhead_SameLane[2][5] > VehAhead_SameLane[0][5]) &&
                             (VehAhead_SameLane[2][5] > VehAhead_SameLane[1][5]))  
                             for(int i=0; i < 7 ; i++) SameLane_abc[2][i] = VehAhead_SameLane[2][i];   

                      else for(int i=0; i < 7 ; i++) SameLane_abc[1][i] = VehAhead_SameLane[2][i];   
                      }

                cout << "SameLane_abc[0]" << SameLane_abc[0][5] << endl;   
                cout << "SameLane_abc[1]" << SameLane_abc[1][5] << endl;   
                cout << "SameLane_abc[2]" << SameLane_abc[2][5] << endl;   
  
                  

  
                int prev_size = previous_path_x.size();


                if (prev_size > 0)
                {
                   car_s = end_path_s;
                }

                bool too_close = false;
  

                double a_vx = SameLane_abc[0][3];
                double a_vy = SameLane_abc[0][4];
                double a_v  = sqrt(a_vx*a_vx+a_vy*a_vy);
                double a_s  = SameLane_abc[0][5];
                a_s += ((double)prev_size*0.02*a_v);
                
                cout <<"diff = "<< a_s-car_s << endl;
 
                if((a_s > car_s) && ((a_s-car_s) < 30))
                       {

                           //ref_vel = 29.5;
                           too_close = true;
                           //close_car_speed = check_speed;
                           if (lane>0)
                           //{
                               lane=0;
                           //}  


                       }
                if(too_close)
                {
                    ref_vel -= 0.1;//diff_speed*0.03;      

 
                }
                else if(ref_vel<49.5)
                {
                    ref_vel += 0.1;//.224; 
                } 




                // Create list of widely space waypoints
                vector<double> ptsx;
                vector<double> ptsy;

                // referencestarting point
                double ref_x = car_x;
                double ref_y = car_y;
                double ref_yaw =  deg2rad(car_yaw);
                
                // if previous size is almost empty, use car as starting reference

                if (prev_size < 2) 
                {  
                    double prev_car_x = car_x - cos(car_yaw);
                    double prev_car_y = car_y - sin(car_yaw);
                      
                    ptsx.push_back(prev_car_x);
                    ptsx.push_back(car_x);
                    
                    ptsy.push_back(prev_car_y);
                    ptsy.push_back(car_y);
                    
                }

                else

                {
                    ref_x = previous_path_x[prev_size-1];
                    ref_y = previous_path_y[prev_size-1]; 

                    double ref_x_prev = previous_path_x[prev_size-2];
                    double ref_y_prev = previous_path_y[prev_size-2];
                    ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);

                    ptsx.push_back(ref_x_prev);
                    ptsx.push_back(ref_x);

                    ptsy.push_back(ref_y_prev);
                    ptsy.push_back(ref_y);
                }


                vector<double> next_mp0 = getXY(car_s+30,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);  
                vector<double> next_mp1 = getXY(car_s+60,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);  
                vector<double> next_mp2 = getXY(car_s+90,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);  

                ptsx.push_back(next_mp0[0]);
                ptsx.push_back(next_mp1[0]);
                ptsx.push_back(next_mp2[0]);       
                
                ptsy.push_back(next_mp0[1]);
                ptsy.push_back(next_mp1[1]);
                ptsy.push_back(next_mp2[1]);

                for(int i = 0; i < ptsx.size();i++)
                {
                    double shift_x = ptsx[i]-ref_x;
                    double shift_y = ptsy[i]-ref_y;

                    ptsx[i] = (shift_x *cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
                    ptsy[i] = (shift_x *sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));

                }   
          	 
                tk::spline s;

                s.set_points(ptsx,ptsy);

                  
                //json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


                for(int i = 0; i < previous_path_x.size(); i++)
                {
                    next_x_vals.push_back(previous_path_x[i]);
                    next_y_vals.push_back(previous_path_y[i]);
                }	       


                double target_x = 30.0; 
                double target_y = s(target_x);
                double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

                double x_add_on = 0; 
             
                for(int i = 1; i <= 50-previous_path_x.size(); i++)
                {
                    double N = (target_dist/(0.02*ref_vel/2.24));
                    double x_point = x_add_on + (target_x)/N;
                    double y_point = s(x_point);

                    x_add_on = x_point;

                    double x_ref = x_point;
                    double y_ref = y_point;

                    x_point = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw));
                    y_point = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw));

                    x_point += ref_x;
                    y_point += ref_y;

                    next_x_vals.push_back(x_point);
                    next_y_vals.push_back(y_point);
                }   
                
                json msgJson;
          	
                // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
                /*  
                //Following for loop added in Walkthrough Video by David Silver
                double dist_inc = 0.3;
                for(int i = 0; i < 50; i++)
                {
                    // David Silver added the follwing 5 lines to make the car follow the lane
                    double next_s = car_s + (i+1)*dist_inc;
                    double next_d = 6;
                    vector<double> xy = getXY(next_s,next_d,map_waypoints_s,map_waypoints_x,map_waypoints_y); 
                    next_x_vals.push_back(xy[0]); 
                    next_y_vals.push_back(xy[1]); 
                    
                    //original starter code, it makes the car move straight in constant speed
                    //next_x_vals.push_back(car_x+(dist_inc*i)*cos(deg2rad(car_yaw)));
                    //next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad(car_yaw)));
                }
                */
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
