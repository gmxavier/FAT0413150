//Tutorial #1
function controlFunction(block)
{
  return -3*block.x;
}

//Tutorial #2
function controlFunction(block)
{
  return -3*block.x -1.5*block.dx;
}

//Tutorial #3
var position_error_integral = 0;
function controlFunction(block)
{
  var delta_t = 0.02; // The simulation time step
  position_error_integral += delta_t * block.x;  
  monitor('position error integral       ', position_error_integral);  
  return -3*block.x -5*block.dx -2*position_error_integral;
}

//Cruise Control Intro
function controlFunction(vehicle){
  return 4 * (vehicle.targetSpeed - vehicle.speed);
}

//Cruise Control 2
function controlFunction(vehicle){
  return 4 * (vehicle.targetSpeed - vehicle.speed);
};

//Cruise Control 2
//This function performs a step test.
//Pick up the data at the browser console.
//Paste it on a Google Sheet and edit it as needed.
//Then use https://pidtuner.com/ to identify the model.
var Tnow = 0;
function controlFunction(vehicle){ 
  var CO = 1;
  if (vehicle.T - Tnow > 1) {
  	console.log(vehicle.T, CO, vehicle.speed)
    Tnow += 1;
  }
  return CO;
};

//Ball on Platform Balance
function controlFunction(ball, piston, hinge, T)
{
  var piston_speed_target = 0;
  if((ball.y - (piston.length - 5.15)) < 1.0 && ball.vy*ball.vy > 0.1)
  {
    piston_speed_target = 0.5 * ball.vy;
  }
  piston_speed_target += 2.0 * (2.0 - piston.length);
  var pistonAcceleration = 40 * (piston_speed_target - piston.speed);
  var hinge_angle_target = 0.25 * ball.vx + 0.2 * ball.x;
  var hinge_speed_target = 10 * (hinge_angle_target - hinge.angle);
  var hingeAcceleration = 40 * (hinge_speed_target - hinge.speed);
  return {pistonAcceleration:pistonAcceleration, hingeAcceleration:hingeAcceleration};
}

//Ball on Platform Bounce
var punt_speed = 3.9;
function controlFunction(ball, piston, hinge, T)
{
  var time_to_impact = -(ball.y + 3.5)/ball.vy;
  var piston_speed_target = -10;
  if(ball.vy < 0 && time_to_impact < 1) piston_speed_target = punt_speed + 0.5 * ball.vy;
  var pistonAcceleration = 40 * (piston_speed_target - piston.speed);
  var hinge_angle_target = 0.08 * ball.vx + 0.06 * ball.x;
  var hinge_speed_target = 10 * (hinge_angle_target - hinge.angle);
  var hingeAcceleration = 40 * (hinge_speed_target - hinge.speed);
  var apogee = ball.y + ball.vy*ball.vy / (2*9.81);
  punt_speed -= 0.001 * apogee;
  monitor('punt_speed   ', punt_speed);
  monitor('apogee       ', ball.y + ball.vy*ball.vy / (2*9.81));
  return {pistonAcceleration:pistonAcceleration, hingeAcceleration:hingeAcceleration};
}

//Rocket Landing
function controlFunction(rocket)
{
  var gimbalAngle1 = 0.01*rocket.x + 0.05*rocket.dx;
  var gimbalAngle2 = 2*rocket.theta + 4*rocket.dtheta;
  var gimbalAngle = gimbalAngle1 + gimbalAngle2;
  var throttle = -0.1*rocket.y - 2*rocket.dy
  return {throttle:throttle, gimbalAngle:gimbalAngle};
}

//Rocket Landing 2
function controlFunction(rocket)
{
  // Horizontal position control
  var SP_x = 0;
  var E_x = SP_x - rocket.x;
  var dE_x = -rocket.dx
  var KP_x = -0.01;
  var KD_x = -0.05;  
  var P_x = KP_x*E_x + KD_x*dE_x;
 
  // Vertical position control
  var SP_y = 10;
  var E_y = SP_y - rocket.y;
  var dE_y = -rocket.dy;
  var KP_y = 0.1;
  var KD_y = 2;
  var P_y = KP_y*E_y + KD_y*dE_y;
  
  // Pitch control
  var SP_theta = 0;
  var E_theta = SP_theta - rocket.theta;
  var dE_theta = -rocket.dtheta;  
  var KP_theta = -2;
  var KD_theta = -4;
  var P_theta = KP_theta*E_theta + KD_theta*dE_theta;
    
  var gimbalAngle = P_x + P_theta;
  var throttle = P_y
  
  return {throttle:throttle, gimbalAngle:gimbalAngle};
}
