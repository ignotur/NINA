#include <iostream>
#include <fstream>

using namespace std;

int main () {
ifstream in ("nevoidN.NE2001.dat");
ofstream out("nevoid.incl.txt");

double value;
double flag;
int n=1;

	do {
		in >> flag;
//		out<<"           voidflag = "<<flag<<endl;
		in >> value;
		out<<"           lv("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           bv("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           dv("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           nev("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           Fv("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           aav("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           bbv("<<n<<") = "<<value<<endl; 
		in >> value;
		out<<"           ccv("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           thvy("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           thvz("<<n<<") = "<<value<<endl;
		in >> value;
		out<<"           edge("<<n<<") = "<<value<<endl; 

		n++;

	} while (!in.eof());

return 0;
}
