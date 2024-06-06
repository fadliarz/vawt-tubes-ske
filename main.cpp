#include <bits/stdc++.h>
#include <cmath>
#define _USE_MATH_DEFINES
#define TO_RAD 0.0174533
#define TO_DEG 1 / 0.0174533
#include <math.h>

using namespace std;

// base parameter
double L1 = 0.070;
double L2 = 0.610;
double L4 = 0.140;
double L3 = 0.55; // initial value
double maximumL3 = 0.68;
double deltaL3 = 0.0001;
double deltaTetha = 1;

int getQuadrant(double tetha) {
  if (0 <= tetha && tetha <= 90) {
    return 1;
  }

  if (90 < tetha && tetha <= 180) {
    return 2;
  }

  if (180 < tetha && tetha <= 270) {
    return 3;
  }

  if (270 < tetha && tetha <= 360) {
    return 4;
  }

  throw 400;
}

double modifyAngleOfAttack(double angleOfAttack) {
  while (angleOfAttack < 0) {
    angleOfAttack += 360;
  }

  while (angleOfAttack > 360) {
    angleOfAttack -= 360;
  }

  return angleOfAttack;
}

pair<double, double> getForce(double tetha, double L3,
                              vector<double> horizontalForce,
                              vector<double> verticalForce,
                              vector<double> horizontalForce2,
                              vector<double> verticalForce2) {

  double M, M2, angleOfAttack, angleOfAttack2;
  double R, T;
  double x, Y, H;
  double B, J;
  double phi;
  double forceQuadrantCoefficient;
  switch (getQuadrant(tetha)) {
  case 1:
  case 2:
    x = sqrt(pow(L1, 2) + pow(L2, 2) - 2 * L1 * L2 * cos(tetha * TO_RAD));
    Y = acos((pow(L4, 2) + pow(x, 2) - pow(L3, 2)) / (2 * L4 * x)) * TO_DEG;
    H = asin(L1 / x * sin(tetha * TO_RAD)) * TO_DEG;
    M = Y + H - 180 + tetha;
    R = asin(L4 / L3 * sin(Y * TO_RAD)) * TO_DEG;
    T = 180 - Y - R;
    angleOfAttack = 360 - M;
    M2 = atan(x / L1 * sin(R * TO_RAD) / sin(T * TO_RAD)) * TO_DEG - T;
    angleOfAttack2 = 360 - M2;
    phi = M;
    forceQuadrantCoefficient = 1;
    break;
  case 3:
  case 4:
    B = 360 - tetha;
    x = sqrt(pow(L1, 2) + pow(L2, 2) - 2 * L1 * L2 * cos(B * TO_RAD));
    Y = acos((pow(L3, 2) + pow(L4, 2) - pow(x, 2)) / (2 * L3 * L4)) * TO_DEG;
    M2 = 90 - Y;
    angleOfAttack2 = -1 * M2;
    J = acos((pow(L4, 2) + pow(x, 2) - pow(L3, 2)) / (2 * L4 * x)) * TO_DEG;
    M = B + J - M2 - 90;
    angleOfAttack = 180 + M;
    forceQuadrantCoefficient = -1;
    break;
  }

  angleOfAttack = modifyAngleOfAttack(angleOfAttack);
  angleOfAttack2 = modifyAngleOfAttack(angleOfAttack2);

  double Fx, Fy;
  int lowerAngleOfAttack = angleOfAttack;
  int upperAngleOfAttack = lowerAngleOfAttack + 1;
  double horizontalForceGradient, verticalForceGradient, Fi, Fj;
  if (upperAngleOfAttack > 359) {
    Fi = horizontalForce[lowerAngleOfAttack];
    Fj = verticalForce[lowerAngleOfAttack];
  } else {
    horizontalForceGradient = (horizontalForce[upperAngleOfAttack] -
                               horizontalForce[lowerAngleOfAttack]) /
                              (upperAngleOfAttack - lowerAngleOfAttack);
    verticalForceGradient = (verticalForce[upperAngleOfAttack] -
                             verticalForce[lowerAngleOfAttack]) /
                            (upperAngleOfAttack - lowerAngleOfAttack);

    Fi = horizontalForceGradient * (angleOfAttack - upperAngleOfAttack) +
         horizontalForce[upperAngleOfAttack];
    Fj = verticalForceGradient * (angleOfAttack - upperAngleOfAttack) +
         verticalForce[upperAngleOfAttack];
  }

  double lowerAngleOfAttack2 = angleOfAttack2;
  double upperAngleOfAttack2 = lowerAngleOfAttack2 + 1;
  double horizontalForceGradient2, verticalForceGradient2, Fi2, Fj2;
  if (upperAngleOfAttack2 > 359) {
    Fi2 = horizontalForce2[lowerAngleOfAttack2];
    Fj2 = verticalForce2[lowerAngleOfAttack2];
  } else {
    horizontalForceGradient2 = (horizontalForce2[upperAngleOfAttack2] -
                                horizontalForce2[lowerAngleOfAttack2]) /
                               (upperAngleOfAttack2 - lowerAngleOfAttack2);
    verticalForceGradient2 = (verticalForce2[upperAngleOfAttack2] -
                              verticalForce2[lowerAngleOfAttack2]) /
                             (upperAngleOfAttack2 - lowerAngleOfAttack2);

    Fi2 = horizontalForceGradient2 * (angleOfAttack2 - upperAngleOfAttack2) +
          horizontalForce2[upperAngleOfAttack2];
    Fj2 = verticalForceGradient2 * (angleOfAttack2 - upperAngleOfAttack2) +
          verticalForce2[upperAngleOfAttack2];
  }

  Fx = forceQuadrantCoefficient *
       ((Fi + Fi2) * cos(phi * TO_RAD) - (Fj + Fj2) * sin(phi * TO_RAD));
  Fy = forceQuadrantCoefficient *
       ((Fi + Fi2) * sin(phi * TO_RAD) + (Fj + Fj2) * cos(phi * TO_RAD));

  pair<double, double> force = make_pair(Fx, Fy);

  pair<double, double> forcePosition =
      make_pair(L2 * cos(tetha * TO_RAD), L2 * sin(tetha * TO_RAD));
  double moment =
      forcePosition.first * force.second - forcePosition.second * force.first;

  return force;
}

int main() {
  ifstream fin;

  // open stream
  fin.open("data.csv");

  if (!fin.is_open())
    cout << "An error occured when opening file" << endl;

  vector<string> row;
  string line, cell;

  vector<double> horizontalForce;
  vector<double> verticalForce;
  while (fin.good()) {
    row.clear();

    getline(fin, line);

    stringstream s(line);

    while (getline(s, cell, ';')) {
      row.push_back(cell);
    }
    horizontalForce.push_back(stod(row[1]));
    verticalForce.push_back(stod(row[2]));
  }

  // open stream
  ifstream fin2;
  fin2.open("F2tambahan.csv");

  if (!fin2.is_open())
    cout << "An error occured when opening file" << endl;

  vector<double> horizontalForce2(360, 0);
  vector<double> verticalForce2(360, 0);
  while (fin2.good()) {
    row.clear();

    getline(fin2, line);

    stringstream s(line);

    while (getline(s, cell, ';')) {
      row.push_back(cell);
    }

    horizontalForce2[stoi(row[1])] = stod(row[2]);
    verticalForce2[stoi(row[1])] = stod(row[3]);
  }

  double currentM;
  double maximumM = INT_MIN;
  double optimumL = -1;

  pair<double, double> force;
  vector<double> vectorOfL3;
  vector<double> vectorOfMoment;
  while (true) {
    currentM = 0;
    for (double tetha = 0; tetha <= 360; tetha = tetha + deltaTetha) {
      force = getForce(tetha, L3, horizontalForce, verticalForce,
                       horizontalForce2, verticalForce2);
      if (!isnan(force.first) && !isnan(force.second)) {
        pair<double, double> forcePosition =
            make_pair(L2 * cos(tetha * TO_RAD), L2 * sin(tetha * TO_RAD));
        currentM += abs(forcePosition.first * force.second -
                        forcePosition.second * force.first);
      } else {
        cout << "Error" << endl;
        throw 400;
      }
    }

    if ((currentM) / 360 > maximumM) {
      maximumM = abs(currentM / 360);
      optimumL = L3;
    }

    vectorOfL3.push_back(L3);
    vectorOfMoment.push_back(currentM / 360);

    L3 += deltaL3;

    if (L3 >= maximumL3) {
      break;
    }
  }

  cout << "Optimum L3: " << optimumL << endl;
  cout << "Optimum M: " << maximumM << endl;

  return 0;
}
