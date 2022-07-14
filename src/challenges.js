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
