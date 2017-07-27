#include "gridfun2d.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "dirs.h"
#include "generatormt.h"

using namespace std;
bool gridfun2d::load(const char* filename) {
  std::ifstream input;
  open_data_file(input, filename);
  if (input.fail()) {
    std::cerr << "Indispensable file '" << filename << "' not found"
              << std::endl;
    return false;
  }

  eRes = pRes = 0;
  eMin = eMax = pMin = pMax = 0.0;

  delete[] table;
  table = NULL;
  delete[] erow;
  erow = NULL;
  delete[] prow;
  prow = NULL;

  input >> eRes >> pRes;
  input >> eMin >> pMin;
  input >> eMax >> pMax;

  if (!input) return false;

  table = new double[eRes * pRes];
  erow = new double[eRes];
  prow = new double[pRes];
  if (!(table && erow && prow)) return false;

  double p = 0, e = 0, val = 0;

  for (int i = 0; i < pRes; i++) {
    input >> p;
    double s = 0;
    for (int j = 0; j < eRes; j++) {
      input >> e >> val;
      table[i * eRes + j] = val;
      s += val;
    }
    prow[i] = s;
  }
  if (input)
    return true;
  else
    return false;
}

double gridfun2d::value(const double p, const double e) const {
  if (p < pMin || p >= pMax || e < eMin || e >= eMax || table == NULL)
    return 0.0;

  double rp = (p - pMin) / (pMax - pMin) * pRes - 0.5;
  double re = (e - eMin) / (eMax - eMin) * eRes - 0.5;

  const int np = floor(rp);
  const int ne = floor(re);

  const double pR = rp - np;
  const double eR = re - ne;

  const double c00 = val(np, ne);          // left lower corner
  const double c10 = val(np + 1, ne);      // right lower corner
  const double c01 = val(np, ne + 1);      // left upper corner
  const double c11 = val(np + 1, ne + 1);  // right upper corner

  double f = (c11 - c10 - c01 + c00) * eR * pR + (c10 - c00) * pR +
             (c01 - c00) * eR + c00;

  return f;
}

double gridfun2d::generateE(const double p) const {
  if (p < pMin || p >= pMax || table == NULL) return 0.0;

  double rp = (p - pMin) / (pMax - pMin) * pRes - 0.5;
  const int np = floor(rp);
  const double pR = rp - np;
  double s = 0;

  for (int i = 0; i < eRes; i++)
    s += (erow[i] = (1 - pR) * val(np, i) + pR * val(np + 1, i));

  //    for(int i=eRes-1;i>0;i--)
  //       erow[i]-=(erow[i]-erow[i-1])/2;
  //    erow[0]/=2;

  double des = frandom() * (s - 0.25 * erow[0] - 0.25 * erow[eRes - 1]);

  double real = 0.25 * erow[0];

  if (des < 0.25 * erow[0])
    return eMin + sqrt(des / erow[0]) * (eMax - eMin) / eRes;

  int i = 0;
  while (real < des and i < eRes) {
    real += (erow[i] + erow[i + 1]) / 2;
    i++;
  }
  // now real>=res
  if (i == eRes)
    return eMax - sqrt((des - real) / erow[eRes - 1]) * (eMax - eMin) / eRes;
  else  // now real>=des
  {
    double a = erow[i - 1];
    double b = erow[i];
    double P = real - des;
    double x;
    if (a == b)
      x = P / a;
    else
      x = (sqrt(b * b + 2 * P * (a - b)) - b) / (a - b);
    return eMin + (eMax - eMin) / eRes * (i - x);
  }
}
