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

using namespace std;

// for convenience
using json = nlohmann::json;
// conversions

#define CONVERTMPH2MS 0.447
#define CONVERTMS2MPH 2.24

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

struct Vehicle{
  double x, y, vx, vy, s, d, yaw, speed, timeInLane;
  int id, lane;
  double getSpeed(){
    return sqrt(this->vx*this->vx + this->vy*this->vy);
  }
  int getLane(){
    return int(this->d/4);
  }
  Vehicle():x(0),y(0),vx(0),vy(0),s(0),d(0),id(-1),lane(-1), yaw(0), speed(0){};
};


class Planner{
public:
  // constructor:
  Planner(){};
  // default destructor:
  ~Planner()= default;
  // intializer:
  void initialize(int lane, int numLanes, double vRef, double vTarget, double dt, const vector<double> &map_x,const vector<double> &map_y,const vector<double> &map_s,const vector<double> &map_dx,const vector<double> &map_dy){
    this->lane = lane;
    this->egoLane = lane;
    this->numLanes = numLanes;
    this->vRef = vRef;
    this->vTarget = vTarget;
    this->dt = dt;
    this->map_waypoints_x = map_x;
    this->map_waypoints_y = map_y;
    this->map_waypoints_s = map_s;
    this->maxAccel = 10;
    this->maxJerk = 10;
  };
  // Processing step:
  void advance(Vehicle &ego, Vehicle &ego_prev, vector<Vehicle> &cars, vector<double> &prevPathX, vector<double> &prevPathY, vector<double> &trajectory_x, vector<double> & trajectory_y){
    // Reset sensor data logicals:
    this->reset();
    // Housekeeping some ego testing, updating lane timers...
    this->updateEgo(ego, ego_prev, !prevPathX.empty());
    // Determine if left/right lanes are clear:
    this->checkLanes(ego, cars);
    // Update FSM state based on relative locations & speeds of surrounding vehicles:
    this->FSMUpdate(ego, cars, prevPathX.size());
    // Generate Trajectory
    this->generateTrajectory(ego, prevPathX, prevPathY, trajectory_x, trajectory_y);
    this->seq++;
  };
private:
  // VARIABLE DEFINITION:
  int lane, numLanes, egoLane;
  // reference velocity
  double vRef, vTarget;
  double maxAccel, maxJerk;
  // Time step:
  double dt;
  // sequence number
  unsigned long seq = 0;
  // Mapping variable:
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;

  // Logicals:
  bool eBrake;
  bool leftClear ;
  bool rightClear;
  bool keepLane;
  // Speed variables:
  double minDistance = 30;
  // positioning variables:
  double forwardHorizon = 30;
  double forwardSteps = 50;

  //PRIVATE FUNCTIONS:
  void reset(){
//    this->is_car_too_close = false;
    this->eBrake = false;
    this->leftClear = true;
    this->rightClear = true;
    this->keepLane = false;
  };

  void updateEgo(Vehicle &ego, Vehicle &ego_prev, bool prevPathExists){
    // get time spent in current lane:
    if (this->seq>0){
      if (this->egoLane == ego.getLane()){
        ego.timeInLane += this->dt;
      }else{
        ego.timeInLane = 0;
      }
    }
    ego.lane = ego.getLane();
    this->egoLane = ego.lane;

    // set car's current position to last path position if it exists:
    if (prevPathExists){
      ego.s = ego_prev.s;
    }
  }

  void checkLanes(Vehicle &ego, vector<Vehicle> &cars){
    for(auto &car: cars) {
      if ((car.d < (2 + 4 * (this->lane - 1) + 2)) && (car.d > (2 + 4 * (this->lane - 1) - 2))) {
        // The car is in the left lane to the vehicle
        if (abs(car.s - ego.s) < this->minDistance) {
          this->leftClear = false;
        }
      }
      else if ((car.d < (2 + 4 * (this->lane + 1) + 2)) && (car.d > (2 + 4 * (this->lane + 1) -2))) {
        // The car is in the right lane to the vehicle
        if (abs(car.s - ego.s) < this->minDistance) {
          this->rightClear = false;
        }
      }
      // make L/R lanes not clear when in the edge lanes
      if (this->lane == 0) {
        this->leftClear = false;
      }
      else if (this->lane == this->numLanes-1) {
        this->rightClear = false;
      }
    }
  };

  void FSMUpdate(Vehicle &ego, vector<Vehicle> &cars, unsigned long prevPathSize){
    // generic variable for speed of car infront...
    double vCarInFront = 0;
    double distanceToCarInFront = 1e6;
    int startLane = this->lane;
    for(auto &car: cars) {
      // For each other car, first check if it is current lane
      if ((car.d < (2 + 4 * this->lane + 2)) && (car.d > (2 + 4 * this->lane - 2))) {
        // set other car's position based on its speed & the timestep:
        car.s += (double)prevPathSize * this->dt * car.getSpeed();
        // Get distance to car in front:
        double check_distance = car.s - ego.s;
        distanceToCarInFront =(check_distance < distanceToCarInFront)? check_distance: distanceToCarInFront;
        // If the car infront is too close, modify current speed & try to change state to lane change:
        if ((car.s > ego.s) && (check_distance < this->minDistance+10)) {
          vCarInFront = car.getSpeed();
          // three choices based on lanes, inside lane, or L/R outside lanes:
          // if in non-edge lane:
          if ((this->lane > 0) && (this->lane < this->numLanes-1)) {
            // Give motion priority to move to left lane
            if (this->leftClear) {
              this->lane--;
            }
            else if (this->rightClear) {
              this->lane++;
            }
              // no lanes are clear, stay in current lane.
            else {
              this->keepLane = true;
            }
          }
            // if in left lane, can only move to right lane:
          else if (this->lane == 0) {
            if (this->rightClear) {
              this->lane++;
            }
            else {
              this->keepLane = true;
            }
          }
            // if in right lane, can only move to left lane:
          else {
            if (this->leftClear) {
              this->lane--;
            }
            else {
              this->keepLane = true;
            }
          }
          // e-brake if unsafe distance!
          this->eBrake = check_distance<(this->minDistance/4.0);
        }
      }
    }

    /*
     * this block of code implements a 'Passing Left Lane' driving ettiquette that will move the car
     * over to the right lane if it is clear to keep left lanes for passing vehicles.
     * Time in lane ensures low jerk transitions when moving over across all lanes..
     * */
    if ((startLane==this->lane) && (this->lane < this->numLanes-1) && (this->rightClear) && (this->vRef > 45) && (ego.timeInLane>3)){
      this->keepLane = false;
      this->lane++;
    }

    if ((this->keepLane) && (ego.speed > vCarInFront) && (distanceToCarInFront < this->minDistance+10)) {

      if(this->eBrake){
        this->vRef -= (this->maxAccel * this->dt + 0.5 * this->maxJerk * this->dt * this->dt)*CONVERTMS2MPH;
      }else{
        // check how close vRef is to vCarInFront:
        if(fabs(this->vRef*CONVERTMPH2MS - vCarInFront) < this->maxAccel * this->dt){
          this->vRef = vCarInFront*CONVERTMS2MPH;
        }else{
          this->vRef -= (0.5 * this->maxAccel * this->dt)*CONVERTMS2MPH;
        }

      }
    }
    else if (this->vRef < this->vTarget  - 0.5 - (this->maxAccel * this->dt - 0.5 * this->maxJerk * this->dt * this->dt)*CONVERTMS2MPH) {
      this->vRef += (this->maxAccel * this->dt + 0.5 * this->maxJerk * this->dt * this->dt)*CONVERTMS2MPH;
    }
  }

  void generateTrajectory(Vehicle &ego,  vector<double> &prevX, vector<double> &prevY, vector<double> &trajectory_x, vector<double> &trajectory_y){
    vector<double> pX;
    vector<double> pY;

    // Current x,y position of ego vehicle:
    double rX = ego.x;
    double rY = ego.y;
    // current yaw in rad:
    double rYaw = deg2rad(ego.yaw);

    // adding previous position of vehicle
    if(prevX.size() < 2) {
      // is previous path has 1 or less points, subtrack points back using current yaw:
      double prev_car_x = ego.x - cos(ego.yaw);
      double prev_car_y = ego.y - sin(ego.yaw);

      pX.emplace_back(prev_car_x);
      pX.emplace_back(ego.x);

      pY.emplace_back(prev_car_y);
      pY.emplace_back(ego.y);
    }
    else {
      // if previous path has 2 or more points, use previous path to fill last points & get the theta:
      rX = prevX[prevX.size() - 1];
      rY = prevY[prevX.size() - 1];

      double ref_x_prev = prevX[prevX.size() - 2];
      double ref_y_prev = prevY[prevX.size() - 2];
      rYaw = atan2(rY - ref_y_prev, rX -ref_x_prev);

      pX.emplace_back(ref_x_prev);
      pX.emplace_back(rX);

      pY.emplace_back(ref_y_prev);
      pY.emplace_back(rY);

    }

    // Add n-points to end of path at 30m increments, thus, the lane change is completed in 30m:
    for (int i=0; i<3; i++){
      vector<double> wp;
      wp = getXY(ego.s + this->forwardHorizon*(i+1), (2 + 4 * this->lane), this->map_waypoints_s, this->map_waypoints_x, this->map_waypoints_y);
      pX.emplace_back(wp[0]);
      pY.emplace_back(wp[1]);
    }

    // convert points into car perspective:
    for(int i = 0; i < pX.size(); i++) {
      // Shifting car's reference angle
      double dx = pX[i] - rX;
      double dy = pY[i] - rY;
      // rotate points to align with camera's reference:
      pX[i] = (dx * cos(0 - rYaw) - dy * sin(0 - rYaw));
      pY[i] = (dx * sin(0 - rYaw) + dy * cos(0 - rYaw));
    }
    // create spline object:
    tk::spline spline;
    // fit spline in car's reference point:
    spline.set_points(pX, pY);

    // add previous trajectory points to current trajectory to smooth:
    for (int i = 0; i < prevX.size(); i++) {
      trajectory_x.emplace_back(prevX[i]);
      trajectory_y.emplace_back(prevY[i]);
    }

    // set forward motion target:
    double tY = spline(this->forwardHorizon);
    double target_dist = sqrt(this->forwardHorizon * this->forwardHorizon + tY * tY);

    double xSum = 0;

    // NOTE: vRef was already set previously, thus, need to take the spline & discretize it so change in distance coordinate with vRef & dt!
    for (int i = 1; i <= this->forwardSteps - prevX.size(); i++) {
      double x = xSum + this->dt*this->vRef*CONVERTMPH2MS; // convert vRef to m/s with minor reduction for cross-lane jerk
      double y = spline(x);

      xSum = x;

      double x_ref = x;
      double y_ref = y;
      // change from vehicle to global reference frame:
      x = x_ref * cos(rYaw) - y_ref * sin(rYaw);
      y = x_ref * sin(rYaw) + y_ref * cos(rYaw);

      x += rX;
      y += rY;

      trajectory_x.emplace_back(x);
      trajectory_y.emplace_back(y);
    }
  }
};

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
    map_waypoints_x.emplace_back(x);
    map_waypoints_y.emplace_back(y);
    map_waypoints_s.emplace_back(s);
    map_waypoints_dx.emplace_back(d_x);
    map_waypoints_dy.emplace_back(d_y);
  }

  // Start in lane 1
  int lane = 1;
  // total number of lanes:
  int numLanes = 3;

  // Reference velocity to target
  double vRef = 0.0; //mph
  // set target speed:
  double vTarget = 50;
  // Time step for simulator:
  double dt = 0.02; //s

  Planner planner;
  planner.initialize(lane, numLanes, vRef, vTarget, dt, map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy);
  h.onMessage([&planner](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,uWS::OpCode opCode) {
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

          // Main car's localization Data
          Vehicle ego;
          ego.x = j[1]["x"];
          ego.y = j[1]["y"];
          ego.s = j[1]["s"];
          ego.d = j[1]["d"];
          ego.yaw = j[1]["yaw"];
          ego.speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // push previous paths into vector<double> types:
          vector<double> _prevPathX, _prevPathY;
          for (int i=0; i < previous_path_x.size(); i++){
            _prevPathX.emplace_back(previous_path_x[i]);
            _prevPathY.emplace_back(previous_path_y[i]);
          }
          // Previous path's end s and d values
          Vehicle ego_prev;
          ego_prev.s = j[1]["end_path_s"];
          ego_prev.d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];
          // put sensor data into vector of vehicle struct:
          vector<Vehicle> vehicles;
          for (int i=0; i < sensor_fusion.size(); i++){
            Vehicle v;
            v.id = sensor_fusion[i][0];
            v.x = sensor_fusion[i][1];
            v.y = sensor_fusion[i][2];
            v.vx = sensor_fusion[i][3];
            v.vy = sensor_fusion[i][4];
            v.s = sensor_fusion[i][5];
            v.d = sensor_fusion[i][6];
            vehicles.emplace_back(v);
          }

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;
          // run update on planner:
          planner.advance(ego, ego_prev, vehicles, _prevPathX, _prevPathY, next_x_vals, next_y_vals);

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
