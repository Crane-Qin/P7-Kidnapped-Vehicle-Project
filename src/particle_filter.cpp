/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#define EPS 0.0001

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  
  default_random_engine gen;
  
  // Creates a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  num_particles = 120;
  
  for (int i = 0; i < num_particles; ++i) {
    
    // TODO: Sample  and from these normal distrubtions like this:
    //   sample_x = dist_x(gen);
    //   where "gen" is the random engine initialized earlier.
    
    Particle particle;
    
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particle.id = i;
    particles.push_back(particle);
  }
  
  is_initialized =true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  
  default_random_engine gen;
  
  // Make distributions for adding noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  // Calculate new state.
  for (int i = 0; i < num_particles; i++) {
    if (fabs(yaw_rate)<EPS){
      particles[i].x +=velocity*delta_t*cos(particles[i].theta);
      particles[i].y +=velocity*delta_t*sin(particles[i].theta);
      //yaw not change
    }
    else{
      particles[i].x +=velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
      particles[i].y +=velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    
    // Add noise to the particles
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.
  for (unsigned int i =0; i< observations.size();i++){
    
    double min_distance = numeric_limits<double>::max();
    int map_id = 0;
    for(unsigned int j =0; j< predicted.size();j++){
      double distance = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
      //double distance = sqrt((observations[i].x-predicted[j].x)*(observations[i].x-predicted[j].x)+(observations[i].y-predicted[j].y)*(observations[i].y-predicted[j].y));
      if (distance<min_distance){
        min_distance = distance;
        map_id = predicted[j].id;
      }
    }
    observations[i].id = map_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks){
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html
  
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  
  for (unsigned int i = 0; i < num_particles; i++) {
    
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    //creat vector to store landmarks in sensor range
    vector<LandmarkObs> validLandmarks;
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      float x_landmark = map_landmarks.landmark_list[j].x_f;
      float y_landmark = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      double distance = sqrt((x-x_landmark)*(x-x_landmark)+(y-y_landmark)*(y-y_landmark));
      if ( distance<= sensor_range ) {
        validLandmarks.push_back(LandmarkObs{ id, x_landmark, y_landmark });
      }
    }
    
    //creat vector to store transformed observations
    vector<LandmarkObs> observationsTOmap;
    for(unsigned int j = 0; j < observations.size(); j++) {
      //Homogenous Transformation : from car coordinate system to map coordinate
      int id =observations[j].id;
      double x_map = x+ cos(theta)*observations[j].x - sin(theta)*observations[j].y ;
      double y_map = y+ sin(theta)*observations[j].x + cos(theta)*observations[j].y ;
      observationsTOmap.push_back(LandmarkObs{ id, x_map, y_map });
    }
    
    //creat vector to store transformed observations
    dataAssociation(validLandmarks, observationsTOmap);
    
    // set weight for each particle
    particles[i].weight = 1.0;
    // calculate weights
    for(unsigned int j = 0; j < observationsTOmap.size(); j++) {
      double x_obs = observationsTOmap[j].x;
      double y_obs = observationsTOmap[j].y;
      
      int id_landmark = observationsTOmap[j].id;
      
      double x_landmark, y_landmark;
      unsigned int k = 0;
      //unsigned int n = validLandmarks.size();
      bool detected = false;
      while( !detected && k < validLandmarks.size() ) {
        if ( validLandmarks[k].id == id_landmark) {
          detected = true;
          x_landmark = validLandmarks[k].x;
          y_landmark = validLandmarks[k].y;
        }
        k++;
      }
      
      double exponent =1.0/2* ((x_obs-x_landmark)*(x_obs-x_landmark)/sig_x/sig_x + (y_obs-y_landmark)*(y_obs-y_landmark)/sig_y/sig_y);
      //accumulate weights
      particles[i].weight *= 1.0/(2*M_PI*sig_x*sig_y)*exp(-exponent);
    }
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  default_random_engine gen;
  
  vector<Particle> new_particles;
  //create vector to store each particle's weight
  vector<double> weights;
  
  for (unsigned int i=0;i<num_particles;i++){
    weights.push_back(particles[i].weight);
  }
  double max_weight= *max_element(weights.begin(), weights.end());
  double beta=0.0;
  
  // get random initial index.
  uniform_int_distribution<int> random_index(0, num_particles - 1);
  int index = random_index(gen);
  
  uniform_real_distribution<double> random_weight(0.0, 1.0);
  
  for (unsigned int i=0;i<num_particles;i++){
    beta += random_weight(gen) * 2.0 * max_weight;
    while(beta>weights[index]){
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }
  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                         const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  
  //clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();
  
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
  return particle;
  
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best)
{
  vector<double> v = best.sense_x;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best)
{
  vector<double> v = best.sense_y;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}



