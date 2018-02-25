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
          	double Five_x = j[1]["x"];  
          	double Five_y = j[1]["y"];
          	double Five_s = j[1]["s"];
          	double Five_d = j[1]["d"];
          	double Five_yaw = j[1]["yaw"];
          	double car_x_speed = j[1]["speed"];

                //double car_a_speed = .0;
                //double Two_sa_s = 0.0;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	
                vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];
               


                //     When 5 is in middle lane 
                //
                //     | 1 | 2 | 3 |  <<< OneTwoThree
                //     | x | 5 | x |   
                //     | 7 | x | 9 |  <<< SevenNine    

                //     When 5 is in right lane
                //
                //     | x | 1 | 2 |  <<< OneTwoThree
                //     | x | x | 5 |   
                //     | x | 7 | x |  <<< SevenNine    

                //     When 5 is in left lane
                //
                //     | 2 | 3 | x |  <<< OneTwoThree
                //     | 5 | x | x |   
                //     | x | 9 | x |  <<< SevenNine    



                // Get OneTwoThree from sensor_fusion
                // a vector of vectors with all vehicles within 400 meters ahead

                vector<vector<double>> OneTwoThree = j[1]["sensor_fusion"];
                     
                int p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {
                    
                    OneTwoThree[i][0] = 0.0;
                    OneTwoThree[i][1] = 0.0;
                    OneTwoThree[i][2] = 0.0;
                    OneTwoThree[i][3] = 0.0;
                    OneTwoThree[i][4] = 0.0;
                    OneTwoThree[i][5] = 0.0;
                    OneTwoThree[i][6] = 0.0;
                    
                    if (((sensor_fusion[i][5] - Five_s) < 400.0) && ((sensor_fusion[i][5] - Five_s) > 0))
                        {
                                         
                        OneTwoThree[p][0] = sensor_fusion[i][0];
                        OneTwoThree[p][1] = sensor_fusion[i][1];
                        OneTwoThree[p][2] = sensor_fusion[i][2];
                        OneTwoThree[p][3] = sensor_fusion[i][3];
                        OneTwoThree[p][4] = sensor_fusion[i][4];
                        OneTwoThree[p][5] = sensor_fusion[i][5];
                        OneTwoThree[p][6] = sensor_fusion[i][6];
                        p++;
                        
                        }
                }


                // Get SevenNine from sensor_fusion
                // a vector of vectors with all vehicles within 400 meters behind



                vector<vector<double>> SevenNine = j[1]["sensor_fusion"];
                     
                p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {
                    
                    SevenNine[i][0] = 0.0;
                    SevenNine[i][1] = 0.0;
                    SevenNine[i][2] = 0.0;
                    SevenNine[i][3] = 0.0;
                    SevenNine[i][4] = 0.0;
                    SevenNine[i][5] = 0.0;
                    SevenNine[i][6] = 0.0;
                    
                    if (((Five_s - sensor_fusion[i][5]) < 400.0) && ((Five_s - sensor_fusion[i][5]) > 0))
                        {
                                         
                        SevenNine[p][0] = sensor_fusion[i][0];
                        SevenNine[p][1] = sensor_fusion[i][1];
                        SevenNine[p][2] = sensor_fusion[i][2];
                        SevenNine[p][3] = sensor_fusion[i][3];
                        SevenNine[p][4] = sensor_fusion[i][4];
                        SevenNine[p][5] = sensor_fusion[i][5];
                        SevenNine[p][6] = sensor_fusion[i][6];
                        p++;
                        
                        }
                }


                // Get One from OneTwoThree
                // a vector of vectors with all vehicles at 400 meters ahead and on left lane

                vector<vector<double>> One = j[1]["sensor_fusion"];
                     
                p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {                    
                    One[i][0] = 0.0;
                    One[i][1] = 0.0;
                    One[i][2] = 0.0;
                    One[i][3] = 0.0;
                    One[i][4] = 0.0;
                    One[i][5] = 0.0;
                    One[i][6] = 0.0;
                    
                    if ((-OneTwoThree[i][6] + Five_d) > 2.0)
                        {                         
                        One[p][0] = OneTwoThree[i][0];
                        One[p][1] = OneTwoThree[i][1];
                        One[p][2] = OneTwoThree[i][2];
                        One[p][3] = OneTwoThree[i][3];
                        One[p][4] = OneTwoThree[i][4];
                        One[p][5] = OneTwoThree[i][5];
                        One[p][6] = OneTwoThree[i][6];
                        p++;
              
                        }
                }


                cout << "One0= "<< One[0][5] << "s  " << One[0][6] << "d " << endl;
                cout << "One1= "<< One[1][5] << "s  " << One[1][6] << "d " << endl;
                cout << "One2= "<< One[2][5] << "s  " << One[2][6] << "d " << endl;
                //cout << "One3= "<< One[3][5] << "s  " << One[3][6] << "d " << endl;
                //cout << "One4= "<< One[4][5] << "s  " << One[4][6] << "d " << endl;
                
                // Get Two from OneTwoThree
                // a vector of vectors with all vehicles at 400 meters ahead and on same lane

                vector<vector<double>> Two = j[1]["sensor_fusion"];
                     
                p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {                    
                    Two[i][0] = 0.0;
                    Two[i][1] = 0.0;
                    Two[i][2] = 0.0;
                    Two[i][3] = 0.0;
                    Two[i][4] = 0.0;
                    Two[i][5] = 0.0;
                    Two[i][6] = 0.0;
                    
                    if ((abs(OneTwoThree[i][6] - Five_d) < 2.0))
                        {                         
                        Two[p][0] = OneTwoThree[i][0];
                        Two[p][1] = OneTwoThree[i][1];
                        Two[p][2] = OneTwoThree[i][2];
                        Two[p][3] = OneTwoThree[i][3];
                        Two[p][4] = OneTwoThree[i][4];
                        Two[p][5] = OneTwoThree[i][5];
                        Two[p][6] = OneTwoThree[i][6];
                        p++;
                         
                        }
                }

               
                // Get Three from OneTwoThree
                // a vector of vectors with all vehicles at 400 meters ahead and on rigth lane

                vector<vector<double>> Three = j[1]["sensor_fusion"];
                     
                p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {                    
                    Three[i][0] = 0.0;
                    Three[i][1] = 0.0;
                    Three[i][2] = 0.0;
                    Three[i][3] = 0.0;
                    Three[i][4] = 0.0;
                    Three[i][5] = 0.0;
                    Three[i][6] = 0.0;
                    
                    if ((OneTwoThree[i][6] - Five_d) > 2.0)
                        {                         
                        Three[p][0] = OneTwoThree[i][0];
                        Three[p][1] = OneTwoThree[i][1];
                        Three[p][2] = OneTwoThree[i][2];
                        Three[p][3] = OneTwoThree[i][3];
                        Three[p][4] = OneTwoThree[i][4];
                        Three[p][5] = OneTwoThree[i][5];
                        Three[p][6] = OneTwoThree[i][6];
                        p++;
              
                        }
                }


                cout << "Three0= "<< Three[0][5] << "s  " << Three[0][6] << "d " << endl;
                cout << "Three1= "<< Three[1][5] << "s  " << Three[1][6] << "d " << endl;
                cout << "Three2= "<< Three[2][5] << "s  " << Three[2][6] << "d " << endl;
                //cout << "Three3= "<< Three[3][5] << "s  " << Three[3][6] << "d " << endl;
                //cout << "Three4= "<< Three[4][5] << "s  " << Three[4][6] << "d " << endl;
  


                // Get Seven from SevenNine
                // a vector of vectors with all vehicles at 400 meters behind and on left lane

                vector<vector<double>> Seven = j[1]["sensor_fusion"];
                     
                p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {                    
                    Seven[i][0] = 0.0;
                    Seven[i][1] = 0.0;
                    Seven[i][2] = 0.0;
                    Seven[i][3] = 0.0;
                    Seven[i][4] = 0.0;
                    Seven[i][5] = 0.0;
                    Seven[i][6] = 0.0; 
 
                    if ((-SevenNine[i][6] + Five_d) > 2.0)
                        {                         
                        Seven[p][0] = SevenNine[i][0];
                        Seven[p][1] = SevenNine[i][1];
                        Seven[p][2] = SevenNine[i][2];
                        Seven[p][3] = SevenNine[i][3];
                        Seven[p][4] = SevenNine[i][4];
                        Seven[p][5] = SevenNine[i][5];
                        Seven[p][6] = SevenNine[i][6];
                        p++;
              
                        }
                }


                //cout << "Seven= "<< Seven[0][5] << "s  " << Seven[0][6] << "d " << endl;
                //cout << "Seven= "<< Seven[1][5] << "s  " << Seven[1][6] << "d " << endl;
                //cout << "Seven= "<< Seven[2][5] << "s  " << Seven[2][6] << "d " << endl;
                //cout << "Seven= "<< Seven[3][5] << "s  " << Seven[3][6] << "d " << endl;
                //cout << "Seven= "<< Seven[4][5] << "s  " << Seven[4][6] << "d " << endl;

                // Get Ninefrom SevenNine
                // a vector of vectors with all vehicles at 400 meters behind and on right lane

                vector<vector<double>> Nine= j[1]["sensor_fusion"];
                     
                p = 0;     
                
                for(int i = 0; i < sensor_fusion.size();i++)
                {                    
                    Nine[i][0] = 0.0;
                    Nine[i][1] = 0.0;
                    Nine[i][2] = 0.0;
                    Nine[i][3] = 0.0;
                    Nine[i][4] = 0.0;
                    Nine[i][5] = 0.0;
                    Nine[i][6] = 0.0;
                    
                    if ((SevenNine[i][6] - Five_d) > 2.0)
                        {                         
                        Nine[p][0] = SevenNine[i][0];
                        Nine[p][1] = SevenNine[i][1];
                        Nine[p][2] = SevenNine[i][2];
                        Nine[p][3] = SevenNine[i][3];
                        Nine[p][4] = SevenNine[i][4];
                        Nine[p][5] = SevenNine[i][5];
                        Nine[p][6] = SevenNine[i][6];
                        p++;
              
                        }
                }


                //cout << "Nine= "<< Nine[0][5] << "s  " << Nine[0][6] << "d " << endl;
                //cout << "Nine= "<< Nine[1][5] << "s  " << Nine[1][6] << "d " << endl;
                //cout << "Nine= "<< Nine[2][5] << "s  " << Nine[2][6] << "d " << endl;
                //cout << "Nine= "<< Nine[3][5] << "s  " << Nine[3][6] << "d " << endl;
                //cout << "Nine= "<< Nine[4][5] << "s  " << Nine[4][6] << "d " << endl;

                // Get One_s_abc from One
                // a vector with telemetry of the first 3 vehicle in order of appereance 
 

                vector<vector<double>> One_s_abc(3,vector<double>(7));

 
                if((One[1][5] == 0) && (One[2][5] == 0))
                    for(int i=0; i < 7 ; i++) One_s_abc[0][i] = One[0][i];   
                else  if(One[2][5] == 0)
                           {    
                           if(One[0][5] < One[1][5])
                               {
                               for(int i=0; i < 7 ; i++) One_s_abc[0][i] = One[0][i]; 
                               for(int i=0; i < 7 ; i++) One_s_abc[1][i] = One[1][i];
                               }   
                           else
                               {  
                               for(int i=0; i < 7 ; i++) One_s_abc[0][i] = One[1][i];
                               for(int i=0; i < 7 ; i++) One_s_abc[1][i] = One[0][i]; 
                               }  
                           }    
                       
                else { 


                      if((One[0][5] < One[1][5])  && 
                         (One[0][5] < One[2][5]))   
                         for(int i=0; i < 7 ; i++) One_s_abc[0][i] = One[0][i];   

                      else if((One[0][5] > One[1][5]) &&
                              (One[0][5] > One[2][5]))  
                              for(int i=0; i < 7 ; i++) One_s_abc[2][i] = One[0][i];   

                      else for(int i=0; i < 7 ; i++) One_s_abc[1][i] = One[0][i];   
               

                      if((One[1][5] < One[0][5])  && 
                         (One[2][5] < One[2][5]))   
                         for(int i=0; i < 7 ; i++) One_s_abc[0][i] = One[1][i];   

                      else if((One[1][5] > One[0][5]) &&
                              (One[1][5] > One[2][5]))  
                              for(int i=0; i < 7 ; i++) One_s_abc[2][i] = One[1][i];   

                      else for(int i=0; i < 7 ; i++) One_s_abc[1][i] = One[1][i];   


                      if((One[2][5] < One[0][5])  && 
                         (One[2][5] < One[1][5]))   
                         for(int i=0; i < 7 ; i++) One_s_abc[0][i] = One[2][i];   

                      else if((One[2][5] > One[0][5]) &&
                             (One[2][5] > One[1][5]))  
                             for(int i=0; i < 7 ; i++) One_s_abc[2][i] = One[2][i];   

                      else for(int i=0; i < 7 ; i++) One_s_abc[1][i] = One[2][i];   
                      }

                cout << "One_s_abc[0]" << One_s_abc[0][5] << endl;   
                cout << "One_s_abc[1]" << One_s_abc[1][5] << endl;   
                cout << "One_s_abc[2]" << One_s_abc[2][5] << endl;   
 
                // Get Two_s_abc from Two
                // a vector with telemetry of the first 3 vehicle in order of appereance 
 

                vector<vector<double>> Two_s_abc(3,vector<double>(7));

 
                if((Two[1][5] == 0) && (Two[2][5] == 0))
                    for(int i=0; i < 7 ; i++) Two_s_abc[0][i] = Two[0][i];   
                else  if(Two[2][5] == 0)
                           {    
                           if(Two[0][5] < Two[1][5])
                               {
                               for(int i=0; i < 7 ; i++) Two_s_abc[0][i] = Two[0][i]; 
                               for(int i=0; i < 7 ; i++) Two_s_abc[1][i] = Two[1][i];
                               }   
                           else
                               {  
                               for(int i=0; i < 7 ; i++) Two_s_abc[0][i] = Two[1][i];
                               for(int i=0; i < 7 ; i++) Two_s_abc[1][i] = Two[0][i]; 
                               }  
                           }    
                       
                else { 


                      if((Two[0][5] < Two[1][5])  && 
                         (Two[0][5] < Two[2][5]))   
                         for(int i=0; i < 7 ; i++) Two_s_abc[0][i] = Two[0][i];   

                      else if((Two[0][5] > Two[1][5]) &&
                              (Two[0][5] > Two[2][5]))  
                              for(int i=0; i < 7 ; i++) Two_s_abc[2][i] = Two[0][i];   

                      else for(int i=0; i < 7 ; i++) Two_s_abc[1][i] = Two[0][i];   
               

                      if((Two[1][5] < Two[0][5])  && 
                         (Two[2][5] < Two[2][5]))   
                         for(int i=0; i < 7 ; i++) Two_s_abc[0][i] = Two[1][i];   

                      else if((Two[1][5] > Two[0][5]) &&
                              (Two[1][5] > Two[2][5]))  
                              for(int i=0; i < 7 ; i++) Two_s_abc[2][i] = Two[1][i];   

                      else for(int i=0; i < 7 ; i++) Two_s_abc[1][i] = Two[1][i];   


                      if((Two[2][5] < Two[0][5])  && 
                         (Two[2][5] < Two[1][5]))   
                         for(int i=0; i < 7 ; i++) Two_s_abc[0][i] = Two[2][i];   

                      else if((Two[2][5] > Two[0][5]) &&
                             (Two[2][5] > Two[1][5]))  
                             for(int i=0; i < 7 ; i++) Two_s_abc[2][i] = Two[2][i];   

                      else for(int i=0; i < 7 ; i++) Two_s_abc[1][i] = Two[2][i];   
                      }

                cout << "Two_s_abc[0]" << Two_s_abc[0][5] << endl;   
                //cout << "Two_s_abc[1]" << Two_s_abc[1][5] << endl;   
                //cout << "Two_s_abc[2]" << Two_s_abc[2][5] << endl;   
  
                // Get Three_s_abc from Three
                // a vector with telemetry of the first 3 vehicle in order of appereance 
 

                vector<vector<double>> Three_s_abc(3,vector<double>(7));

 
                if((Three[1][5] == 0) && (Three[2][5] == 0))
                    for(int i=0; i < 7 ; i++) Three_s_abc[0][i] = Three[0][i];   
                else  if(Three[2][5] == 0)
                           {    
                           if(Three[0][5] < Three[1][5])
                               {
                               for(int i=0; i < 7 ; i++) Three_s_abc[0][i] = Three[0][i]; 
                               for(int i=0; i < 7 ; i++) Three_s_abc[1][i] = Three[1][i];
                               }   
                           else
                               {  
                               for(int i=0; i < 7 ; i++) Three_s_abc[0][i] = Three[1][i];
                               for(int i=0; i < 7 ; i++) Three_s_abc[1][i] = Three[0][i]; 
                               }  
                           }    
                       
                else { 


                      if((Three[0][5] < Three[1][5])  && 
                         (Three[0][5] < Three[2][5]))   
                         for(int i=0; i < 7 ; i++) Three_s_abc[0][i] = Three[0][i];   

                      else if((Three[0][5] > Three[1][5]) &&
                              (Three[0][5] > Three[2][5]))  
                              for(int i=0; i < 7 ; i++) Three_s_abc[2][i] = Three[0][i];   

                      else for(int i=0; i < 7 ; i++) Three_s_abc[1][i] = Three[0][i];   
               

                      if((Three[1][5] < Three[0][5])  && 
                         (Three[2][5] < Three[2][5]))   
                         for(int i=0; i < 7 ; i++) Three_s_abc[0][i] = Three[1][i];   

                      else if((Three[1][5] > Three[0][5]) &&
                              (Three[1][5] > Three[2][5]))  
                              for(int i=0; i < 7 ; i++) Three_s_abc[2][i] = Three[1][i];   

                      else for(int i=0; i < 7 ; i++) Three_s_abc[1][i] = Three[1][i];   


                      if((Three[2][5] < Three[0][5])  && 
                         (Three[2][5] < Three[1][5]))   
                         for(int i=0; i < 7 ; i++) Three_s_abc[0][i] = Three[2][i];   

                      else if((Three[2][5] > Three[0][5]) &&
                             (Three[2][5] > Three[1][5]))  
                             for(int i=0; i < 7 ; i++) Three_s_abc[2][i] = Three[2][i];   

                      else for(int i=0; i < 7 ; i++) Three_s_abc[1][i] = Three[2][i];   
                      }

                cout << "Three_s_abc[0]" << Three_s_abc[0][5] << endl;   
                cout << "Three_s_abc[1]" << Three_s_abc[1][5] << endl;   
                cout << "Three_s_abc[2]" << Three_s_abc[2][5] << endl;   
                   
                int prev_size = previous_path_x.size();


                if (prev_size > 0)
                {
                   Five_s = end_path_s;
                }

                
  
                double Two_sa_vx = Two_s_abc[0][3];
                double Two_sa_vy = Two_s_abc[0][4];
                double Two_sa_v  = sqrt(Two_sa_vx*Two_sa_vx+Two_sa_vy*Two_sa_vy);
                double Two_sa_s  = Two_s_abc[0][5];
                
                // gap that a vehicle following car a should keep  
                
                double Two_sa_gap = Two_sa_v * 0.4697 * 2.0 ; // 2 seconds in distance
                double Two_sa_rear = Two_sa_s - Two_sa_gap; // s where should Five_s should be   
                double FiveOverTwo_saRear = Two_sa_rear - Five_s; // overlap of Five over Two_sa_rear

                // set ref_vel so that vehicle oscillates 15 meters around gap          
       
                if (((Two_sa_rear == 0) && (ref_vel <45)) || ((FiveOverTwo_saRear > 15)&&(ref_vel<45))) 
                    ref_vel += 0.1;
                else if (FiveOverTwo_saRear < 15) ref_vel += (0.1 * FiveOverTwo_saRear/15); 
                

                // Determine which lane is best


                // If 


                int One_val = 0;
                int Two_val = 0;
                int Three_val = 0;
               
  
                if (One_s_abc[0][5] == 0) One_val = 1;
                if (Two_s_abc[0][5] == 0) Two_val = 1;
                if (Three_s_abc[0][5] == 0) Three_val = 1;
                
                cout << "One_val = " << One_val << endl;
                cout << "Two_val = " << Two_val << endl;
                cout << "Three_val = " << Three_val << endl;                   



                // if laneTwo is better than LaneOne and LaneThree, do not change

                // if current lane=1 and LaneOne and LaneThree are equally better than Lane Two, lane=0 



                // Do not go into into opossite side of the road  
                // if current lane=0, and LaneOne is better, do not change lane

                // Do not go out of the road
                // if current lane=2, and laneThree is better, do not change lane


                // What makes 5 change to lane1 if current Lane is 0?
                //    what should make me stay? 

                if ((lane==0) && (Three_val > Two_val)) lane = 1;

                  

                // What makes 5 change to lane0 or lane2 if curent Lane is 1?                

                if ((lane==1) && (One_val > Two_val) && (One_val>Three_val)) lane = 0;
                else if ((lane==1) && (One_val > Two_val) && (Three_val > Two_val) && (One_val == Three_val)) lane = 0;
                else if ((lane==1) && (Three_val>Two_val) && (Three_val > One_val)) lane = 2;                            

                // What makes 5 to change to lane1 if current Lane is 2?                   

                if ((lane == 2) && (One_val > Two_val)) lane = 1;
           

                

                cout << "LANE " << lane << endl; 




                // Create list of widely space waypoints
                vector<double> ptsx;
                vector<double> ptsy;

                // referencestarting point
                double ref_x = Five_x;
                double ref_y = Five_y;
                double ref_yaw =  deg2rad(Five_yaw);
                
                // if previous size is almost empty, use car as starting reference

                if (prev_size < 2) 
                {  
                    double prev_Five_x = Five_x - cos(Five_yaw);
                    double prev_Five_y = Five_y - sin(Five_yaw);
                      
                    ptsx.push_back(prev_Five_x);
                    ptsx.push_back(Five_x);
                    
                    ptsy.push_back(prev_Five_y);
                    ptsy.push_back(Five_y);
                    
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


                vector<double> next_mp0 = getXY(Five_s+30,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);  
                vector<double> next_mp1 = getXY(Five_s+60,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);  
                vector<double> next_mp2 = getXY(Five_s+90,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);  

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
