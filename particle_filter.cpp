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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
        num_particles = 10;

       	default_random_engine gen;
	
	
	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(x, std[0]);
	
	// TODO: Create normal distributions for y and theta.
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

       	for (int i = 0; i < num_particles; i++) {
		
             Particle p;
             p.id = i;
             p.x = dist_x(gen);
             p.y = dist_y(gen);
             p.theta = dist_theta(gen);
             p.weight = 1;
             particles.push_back(p);
       }
       is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
       default_random_engine gen;
       double xf, yf, thetaf, sin_value, cos_value;

       xf = 0; yf = 0; thetaf = 0;

       for (int i = 0; i < num_particles; i++) {
            if (fabs(yaw_rate) > 0.00001) {
                thetaf = particles[i].theta + (yaw_rate * delta_t);
                sin_value = sin(thetaf) - sin(particles[i].theta);
                cos_value = cos(thetaf) - cos(particles[i].theta);

                xf = particles[i].x + ((velocity/(yaw_rate)) *  sin_value);
                yf = particles[i].y + ((velocity/(yaw_rate)) *  cos_value);
            }
            else {
                thetaf = particles[i].theta;
                xf = particles[i].x + (velocity * delta_t * cos(thetaf));
                yf = particles[i].y + ((velocity * delta_t * sin(thetaf)));
            }
	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(xf, std_pos[0]);
	
	// TODO: Create normal distributions for y and theta.
	normal_distribution<double> dist_y(yf, std_pos[1]);
	normal_distribution<double> dist_theta(thetaf, std_pos[2]);

        
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen) ;

        // cout << "Prediction " << particles[i].x << ", " << particles[i].y << ", " << particles[i].theta << " )" << endl;
       }
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

        LandmarkObs near_landm;
        double x_map, y_map, min_dist, dist_calc, updt_weights;
        double gauss_norm= (1/(2.0 * M_PI * std_landmark[0] * std_landmark[1]));
        double stdl_0 = pow(std_landmark[0],2);
        double stdl_1 = pow(std_landmark[1],2);
        double wgt_prod; 


        // clear weights vector
        weights.clear();

        for (int i = 0; i < int(particles.size()); i++) {
           
            updt_weights = 1.0;
            for (int j = 0; j < int(observations.size()); ++j){
                //# transform to map x coordinate
                x_map = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);

                //# transform to map y coordinate
                y_map= particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);

                //  
                min_dist = sensor_range;


            // Finding best landmark

            for (int k = 0; k < int(map_landmarks.landmark_list.size()); ++k) {

                    if  (((fabs(map_landmarks.landmark_list[k].x_f - x_map)) < sensor_range) && (fabs(map_landmarks.landmark_list[k].y_f - y_map) < sensor_range)) {
                    dist_calc = dist(x_map, y_map, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
//                    if ((dist_calc <= min_dist) || (k == 0)) {
                    if ((dist_calc < min_dist)) {
		       near_landm.y = map_landmarks.landmark_list[k].y_f;
		       near_landm.x = map_landmarks.landmark_list[k].x_f;
		       near_landm.id = map_landmarks.landmark_list[k].id_i;
                       min_dist = dist_calc;
                       cout << " Inside Update near.x " << near_landm.x  << " near.y " << near_landm.y << " k " << k << " min_dist " << min_dist << endl;
                    }
                }
            }

           cout << " Update near.x " << near_landm.x  << " near.y " << near_landm.y << " k " << " min_dist " << min_dist << endl;
           // calculating error and 
           // # calculate normalization term

           double x_obs_mu = (x_map - near_landm.x);
           double y_obs_mu = (y_map - near_landm.y);
          // # calculate exponent
          double exponent= ((x_obs_mu * x_obs_mu)/(2.0 * stdl_0)) + ((y_obs_mu * y_obs_mu)/(2 * stdl_1));

          // # calculate weight using normalization terms and exponent
          wgt_prod= gauss_norm * exp(-exponent); 
          cout << " Update x_map " << x_map  << " y_map " << y_map << " j " << j << " i " << i << endl;
          cout << " Update near.x " << near_landm.x  << " near.y " << near_landm.y << endl;
          cout << " Update _wgt  wgt_prod " << wgt_prod  << " )" << endl;
          updt_weights *= wgt_prod; 
          if (updt_weights == 0.0) {

              cout << " Update _wgt == 0 x_obs_mu " << x_obs_mu  << " )" << endl;
              cout << " Update _wgt == 0 y_obs_mu " << y_obs_mu <<  " )" << endl;
              cout << " Update _wgt == 0 exponent " << exponent <<  " )" << endl;
              cout << " Update _wgt == 0 gauss_nr " << gauss_norm << " )" << endl;

          } 
        }
        particles[i].weight = updt_weights;
        weights.push_back(particles[i].weight);
        // cout << " Update " << particles[i].weight << " )" << endl;
        
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution


        std::vector<Particle> particles_new;

        int index;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::discrete_distribution<int> d(weights.begin(), weights.end());

        for(int n=0; n<num_particles; ++n) {
            // Particle particle_res = particles[d(gen)];
            // particles_new.push_back(particle_res);
            index = d(gen);
            particles_new.push_back(particles[index]);
        }


        particles = particles_new;
}

       

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
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
