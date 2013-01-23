#include <iostream>
#include <fstream>

using namespace std;

int main () {
ifstream in ("neclumpN.NE2001.dat");
ofstream out("neclump.incl.txt");

double value;
double flag;
int n=1;

	do {
		in >> flag;
		in >> value;
		out<<"           lc("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           bc("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           nec("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           Fc("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           dc("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           rc("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           edge("<<n<<") = "<<value<<endl; 

		n++;

	} while (!in.eof());

return 0;
}
