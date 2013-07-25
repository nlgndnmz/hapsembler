/*****************************************************************************
    $Author: $
    $Date: $

	Part of Hapsembler package. See the README file for more information.
    Copyright (C) 2011,  Nilgun Donmez <nild@cs.toronto.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 ****************************************************************************/

#include <cmath>

#include "Naive.h"

Naive::Naive(int max_quality, int over, int max_steps, double prior_of_B)
{
	max_qual = max_quality;
	over -= 1;

	num_steps = max_steps;
	prior = new double[num_steps];		// denotes the belief that the two reads are coming from the same location
	double step = 1.0 / num_steps;

	for(int i=0; i<num_steps; i++)
		prior[i] = step * i;

	p_correct = new double[max_qual];
	p_error = new double[max_qual];

	log_correct = new double[max_qual];
	log_error = new double[max_qual];

	log_same = new double*[max_qual];
	log_diff = new double*[max_qual];

	prior_same = new double*[max_qual];
	prior_diff = new double*[max_qual];
	not_prior = new double*[max_qual];

	for(int i=0; i<max_qual; i++)
	{
		int s = (i<1) ? 1 : i;	// if quality score is 0, pretend it is 1

		p_correct[i] = (1 - pow(0.1, s/10.0));
		p_error[i] = pow(0.1, s/10.0)/over;

		log_correct[i] = 0 - log10( 1 - pow(0.1, s/10.0) );
		log_error[i] = 0 - log10( pow(0.1, s/10.0)/over );

		log_same[i] = new double[num_steps];
		log_diff[i] = new double[num_steps];

		prior_same[i] = new double[max_qual];
		prior_diff[i] = new double[max_qual];
		not_prior[i] = new double[max_qual];

		for(int j=0; j<num_steps; j++)
		{
			log_same[i][j] = log_correct[i];
			log_diff[i][j] = 0 - log10( p_error[i] * prior[j] + p_correct[i] * (1 - prior[j]) );
		}
	}

	for(int i=0; i<max_qual; i++)
	{
		for(int k=0; k<max_qual; k++)
		{
			prior_same[i][k] = log10( p_correct[i]*p_correct[k] + over*p_error[i]*p_error[k] );
			prior_diff[i][k] = log10( p_correct[i]*p_error[k] + p_error[i]*p_correct[k] + (over-1)*p_error[i]*p_error[k] );
			not_prior[i][k] = 0 - log_correct[i] - log_correct[k];
		}
	}

	set_prior(prior_of_B);
}

Naive::Naive(const Naive & nv)
{
	this->max_qual = nv.max_qual;
	this->num_steps = nv.num_steps;

	prior = new double[num_steps];
	double step = 1.0 / num_steps;

	for(int i=0; i<num_steps; i++)
		prior[i] = step * i;

	p_correct = new double[max_qual];
	p_error = new double[max_qual];

	log_correct = new double[max_qual];
	log_error = new double[max_qual];

	log_same = new double*[max_qual];
	log_diff = new double*[max_qual];

	prior_same = new double*[max_qual];
	prior_diff = new double*[max_qual];
	not_prior = new double*[max_qual];

	for(int i=0; i<max_qual; i++)
	{
		p_correct[i] = nv.p_correct[i];
		p_error[i] = nv.p_error[i];

		log_correct[i] = nv.log_correct[i];
		log_error[i] = nv.log_error[i];

		log_same[i] = new double[num_steps];
		log_diff[i] = new double[num_steps];

		prior_same[i] = new double[max_qual];
		prior_diff[i] = new double[max_qual];
		not_prior[i] = new double[max_qual];

		for(int j=0; j<num_steps; j++)
		{
			log_same[i][j] = nv.log_correct[i];
			log_diff[i][j] = nv.log_diff[i][j];
		}
	}

	for(int i=0; i<max_qual; i++)
	{
		for(int k=0; k<max_qual; k++)
		{
			prior_same[i][k] = nv.prior_same[i][k];
			prior_diff[i][k] = nv.prior_diff[i][k];
			not_prior[i][k] = nv.not_prior[i][k];
		}
	}

	probB = nv.probB;
	priorB = nv.priorB;
	priorNotB = nv.priorNotB;
}

void Naive::set_prior(double prior_of_B)
{
	probB = int(prior_of_B * num_steps);
	priorB = log10( prior_of_B );
	priorNotB = log10(1.0 - prior_of_B);
}

Naive::~Naive()
{
	delete [] p_correct;
	delete [] p_error;
	delete [] log_correct;
	delete [] log_error;
	delete [] prior;

	for(int i=0; i<max_qual; i++)
	{
		delete [] log_same[i];
		delete [] log_diff[i];
		delete [] prior_same[i];
		delete [] prior_diff[i];
		delete [] not_prior[i];
	}

	delete [] log_same;
	delete [] log_diff;
	delete [] prior_same;
	delete [] prior_diff;
	delete [] not_prior;
}

