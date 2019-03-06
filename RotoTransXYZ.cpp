// Lara Querciagrossa

#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
#include <string.h>
#define THIRD 0.333333333333333
#define ROOTTHREE 1.73205080756888       

using namespace std;


int main(int argc, char *argv[]) 
{
  const double PI = 3.14159265358979;

  //Initializing output files
  ofstream fileout("rototrans.xyz", ios::out);
  ofstream fout("origcom.xyz", ios::out);

  std::cout << " Welcome to RotoTransXYZ code " << endl;
  std::cout << "  " << endl;
  std::cout << " Usage:  " << endl;
  std::cout << "         ./ RotoTransXYZ <file.xyz> <options>  " << endl;
  std::cout << "  " << endl;
  std::cout << "         options: -t x y z OR -r axis degrees" << endl;
  std::cout << "  " << endl;
  std::cout << " Output files: " << endl;
  std::cout << "               origcom.xyz original molecule centered in (0,0,0)" << endl;
  std::cout << "               rototrans.xyz rototed or translated molecule" << endl;
  std::cout << "  " << endl;
  std::cout << "  " << endl;

  /***************************************  
             Read input file
  ***************************************/
  ifstream filein(argv[1]);

  //Check if input file exists
  if ( !filein.is_open() )
    {
      std::cout << "ERROR: Could not open file " << argv[1] << std::endl;
      return 0;
    }

  //Reading header and start writing output files
  int totatom;
  string null;
  filein >> totatom;
  filein >> null;
  std::cout << " " << std::endl;
  std::cout << "Number of atoms: " << totatom << std::endl;
  std::cout << " " << std::endl;

  string *specie = new string[totatom];
  double *x = new double[totatom];
  double *y = new double[totatom];
  double *z = new double[totatom];

  //Center of mass (c.o.m.) initialization
  double comx, comy, comz;
  comx = 0.0;
  comy = 0.0;
  comz = 0.0;

  //Reading input file (atom names and positions)
  for (int i=0; i<totatom; i++){
    filein >> specie[i] >>  x[i] >> y[i] >> z[i];
    //std::cout << "  " << x[i] << " " << y[i] << " " << z[i] << std::endl;
    comx = comx + x[i];
    comy = comy + y[i];
    comz = comz + z[i];
  }

  //Computing c.o.m.
  comx = comx/totatom;
  comy = comy/totatom;
  comz = comz/totatom;

  //Writing centered configuration to "origcom.xyz" file
  fout << totatom << std::endl;
  fout << "# " << std::endl;      
  for (int i=0; i<totatom; i++)
    fout << specie[i] << " " << x[i]-comx << " " << y[i]-comy << " " << z[i]-comz << std::endl;

  //Write rotated or translated configuration
  double *my_x = new double[totatom];
  double *my_y = new double[totatom];
  double *my_z = new double[totatom];

  //Choose between translation and rotation
  string transl = "-t";
  string rotat = "-r";

  //Translate
  if (argv[2] == transl)
    {
      double x_tr = atof(argv[3]);
      double y_tr = atof(argv[4]);
      double z_tr = atof(argv[5]);
      std::cout << "Translation of " << x_tr << " " << y_tr << " " << z_tr << std::endl;
      for (int i=0; i<totatom; i++)
	{
	  my_x[i] = x[i] + x_tr;
	  my_y[i] = y[i] + y_tr;
	  my_z[i] = z[i] + z_tr;
	}
    }
  //Rotate
  else if (argv[2] == rotat)
    {
      string rotax = argv[3];
      double deg = atof(argv[4]);

      std::cout << "Rotation of " << deg << "° around " << rotax << " axis" << std::endl;

      string xax = "x";
      string yax = "y";
      string zax = "z";

      double angle = deg*PI/180.0;

      //Creating rotation matrix
      //  a11 a12 a13
      //  a21 a22 a23
      //  a31 a32 a33

      double a11, a12, a13;
      double a21, a22, a23;
      double a31, a32, a33;

      if (rotax == xax)
	{
	  //Rotate around x axis
	  
	  a11 = 1.0;
	  a12 = 0.0;
	  a13 = 0.0;

	  a21 = 0.0;
	  a22 = cos(angle);
	  a23 = sin(angle);

	  a31 = 0.0;
	  a32 = -sin(angle);
	  a33 = cos(angle);
	}
      else if (rotax == yax)
	{
	  //Rotate around y axis

	  a11 = cos(angle);
	  a12 = 0.0;	  
	  a13 = -sin(angle);

	  a21 = 0.0;
	  a22 = 1.0;
	  a23 = 0.0;

	  a31 = sin(angle);
	  a32 = 0.0;
	  a33 = cos(angle);
	}
      else if (rotax == zax)
	{
	  //Rotate around z axis

	  a11 = cos(angle);
	  a12 = sin(angle);
	  a13 = 0.0;	  

	  a21 = -sin(angle);
	  a22 = cos(angle);
	  a23 = 0.0;	   

	  a31 = 0.0;
	  a32 = 0.0;
	  a33 = 1.0;
	}
      else
	{
	  std::cout << "ERROR! Invalid axis of rotation" << std::endl;
	  return 0;
	}

      double mycomx, mycomy, mycomz;
      mycomx = 0.0;
      mycomy = 0.0;
      mycomz = 0.0;

      //Rotate positions =
      //rotation matrix * position vector
      for (int i=0; i<totatom; i++)
	{
	  my_x[i] = a11*x[i] + a12*y[i] + a13*z[i];
	  my_y[i] = a21*x[i] + a22*y[i] + a23*z[i];
	  my_z[i] = a31*x[i] + a32*y[i] + a33*z[i];

	  mycomx = mycomx + my_x[i];
	  mycomy = mycomy + my_y[i];
	  mycomz = mycomz + my_z[i];
	}

      mycomx = mycomx/totatom;
      mycomy = mycomy/totatom;
      mycomz = mycomz/totatom;

      //Centering rotated/translated configuration
      for (int i=0; i<totatom; i++)
	{
	  my_x[i] = my_x[i] - mycomx;
	  my_y[i] = my_y[i] - mycomy;
	  my_z[i] = my_z[i] - mycomz;
	}
    }
  //Write rotated/translated configuration
  fileout << totatom << std::endl;
  fileout << "# " << std::endl;      
  for (int i=0; i<totatom; i++)
    fileout << setprecision(10) << specie[i] << " " << my_x[i] << " " << my_y[i] << " " << my_z[i] << std::endl;
 
  std::cout << " " << std::endl;
  std::cout << "DONE!" << std::endl;
  std::cout << " " << std::endl;
    
  filein.close();


  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "            /\\____/\\    " << std::endl;
  std::cout << "           ( 0    = )   " << std::endl;
  std::cout << "            \\ >.<  /    " << std::endl;
  std::cout << "       ,-,  /      \\    " << std::endl;
  std::cout << "       : : |       |    " << std::endl;
  std::cout << "       : :  \\      /    " << std::endl;
  std::cout << "        \\::::|_||_|     " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  //  

  return 0;
}


