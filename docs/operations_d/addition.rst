********************************
Addition of coordinate functions
********************************

Description
===========

Let :math:`F=(f_1,\ldots,f_{m_1}) \in \funct{F}_{n,m_1}`, :math:`G=(g_1,\ldots,g_{m_2}) \in \funct{F}_{n,m_2}` and the function conformed by adding the coordinate functions :math:`(F,G) =(f_1,\dots,f_{m_1},g_1,\ldots,g_{m_2}) \in \funct{F}_{n,m_1+m_2}`. Let :math:`\vec{v} \in \gf{V_{m_1+m_2}}`, :math:`\vec{v_F} \in \gf{V_{m_1}}` and :math:`\vec{v_G} \in \gf{V_{m_2}}` so that :math:`\vec{v} = (\vec{v_F}, \vec{v_G})`. 

.. image:: /images/AddImage.png
   :width: 750 px
   :align: center

Library
=======

This construction can be obtained with the following method:

.. code-block:: c

	void addimage(VBF& X, VBF& F, VBF& G)  

Example
-------

The following program provides the Truth Tables of the different intermediate constructions that allow to obtain CLEFIA :math:`S_0` :math:`8 \times 8` S-box from the Truth Tables of the four 4-bit S-boxes :math:`SS_0,SS_1,SS_2` and :math:`SS_3` in which it is constructed and the Truth Table of the multiplication operation in *0x2* performed in :math:`\gf{GF(2^4)}` defined by the
primitive polynomial :math:`x^4 + x + 1`.

.. code-block:: c

	#include <iostream>
	#include <fstream>
	#include "VBF.h"

	int main(int argc, char *argv[])
	{
	   using namespace VBFNS;

	   VBF F,G,T20,T21,U0,U1,Y0,Y1,Y;
	   NTL::mat_GF2 TSS0, TSS1, TSS2, TSS3, Tmul2;
	   NTL::mat_GF2 T2t0, T2t1, Tu0, Tu1, Ty0, Ty1, Ty;

	   ifstream inputSS0("SS0.tt");
	   if(!inputSS0) {
	      cerr << "Error opening " << "SS0.tt" << endl;
	      return 0;
	   }
	   inputSS0 >> TSS0;
	   inputSS0.close();

	   ifstream inputSS1("SS1.tt");
	   if(!inputSS1) {
	      cerr << "Error opening " << "SS1.tt" << endl;
	      return 0;
	   }
	   inputSS1 >> TSS1;
	   inputSS1.close();

	   ifstream inputSS2("SS2.tt");
	   if(!inputSS2) {
	      cerr << "Error opening " << "SS2.tt" << endl;
	      return 0;
	   }
	   inputSS2 >> TSS2;
	   inputSS2.close();

	   ifstream inputSS3("SS3.tt");
	   if(!inputSS3) {
	      cerr << "Error opening " << "SS3.tt" << endl;
	      return 0;
	   }
	   inputSS3 >> TSS3;
	   inputSS3.close();

	   ifstream inputmul2("Mul2.tt");
	   if(!inputmul2) {
	      cerr << "Error opening " << "Mul2.tt" << endl;
	      return 0;
	   }
	   inputmul2 >> Tmul2;
	   inputmul2.close();

	   cout << "t0=" << endl;
	   cout << TSS0 << endl << endl;
	   cout << "t1=" << endl;
	   cout << TSS1 << endl << endl;
	   F.puttt(TSS1);
	   G.puttt(Tmul2);
	   Comp(T21,F,G);
	   T2t1 = TT(T21);
	   cout << "0x2.t1=" << endl;
	   cout << T2t1 << endl;
	   F.kill();
	   G.kill();
	   F.puttt(TSS0);
	   G.puttt(Tmul2);
	   Comp(T20,F,G);
	   T2t0 = TT(T20);
	   cout << "0x2.t0=" << endl;
	   cout << T2t0 << endl;
	   cout << "u0=t0+0x2.t1=" << endl;
	   F.kill();
	   F.puttt(TSS0);
	   directsum(U0,F,T21);
	   Tu0 = TT(U0);
	   cout << Tu0 << endl;
	   G.kill();
	   cout << "u1=0x2.t0+t1=" << endl;
	   G.puttt(TSS1);
	   directsum(U1,T20,G);
	   Tu1 = TT(U1);
	   cout << Tu1 << endl;
	   G.kill();
	   cout << "y0=SS2(u0)=" << endl;
	   G.puttt(TSS2);
	   Comp(Y0,U0,G);
	   Ty0 = TT(Y0);
	   cout << Ty0 << endl;
	   G.kill();
	   cout << "y1=SS3(u1)=" << endl;
	   G.puttt(TSS3);
	   Comp(Y1,U1,G);
	   Ty1 = TT(Y1);
	   cout << Ty1 << endl;
	   addimage(Y,Y0,Y1);
	   Ty = TT(Y);
	   cout << "y=(y0,y1)=" << endl;
	   cout << Ty << endl;

	   return 0;
	}

The output of this program is described in CLEFIA section in "Analysis of CRYPTEC project cryptographic algorithms".

Note that the output of :math:`S_0` S-box :math:`\vec{y} \in \funct{F}_{8,8}` is defined by the addition of coordinate functions of both :math:`\vec{y_0} \in \funct{F}_{8,4}` and :math:`\vec{y_1} \in \funct{F}_{8,4}`. 