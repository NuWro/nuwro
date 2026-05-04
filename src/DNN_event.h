#ifndef DNN_EVENT_H
#define DNN_EVENT_H

#include "../data/F1F209/F1F209Wrapper.hh"
#include "Utilities.h"
#include "e_el_event.h"
#include "e_spp_event.h"
#include "event1.h"
#include "kinematics.h"
#include "nucleus.h"
#include "particle.h"
#include "pdg.h"
#include "sfevent.h"

#include <cmath>
#include <iostream>
#include <stdio.h>

void DNN_event(params &p, event &e, nucleus &t, F1F209Wrapper &BostedFit,
               bool nc);

void Generate_outgoing_nucleons(params &p, event &e, nucleus &t);

#endif
