/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2008-03-27 19:41:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SelfCenteringRub.cpp,v $

// Written: JAE
// Created: Oct 2007
// Revision: A
//
// Description: This file contains the class implementation for
// SelfCenteringRub.

#include <SelfCenteringRub.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>
#include <elementAPI.h>

void* OPS_SelfCenteringRub()
{
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 5) {
		opserr << "WARNING: Insufficient arguements\n";
		opserr << "Want: uniaxialMaterial SelfCentering tag? k1? k2? k3? ";
		opserr << "ActF? beta? <SlipDef? BearDef? rBear?>" << endln;
		return 0;
	}

	int tag;
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid tag\n";
		return 0;
	}

	double data[7] = { 0,0,0,0,0,0,0 };
	numdata = OPS_GetNumRemainingInputArgs();
	if (numdata > 7) {
		numdata = 7;
	}
	if (OPS_GetDoubleInput(&numdata, data)) {
		opserr << "WARNING invalid double inputs\n";
		return 0;
	}

	UniaxialMaterial* mat = new SelfCenteringRub(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]);
	if (mat == 0) {
		opserr << "WARNING: failed to create SelfcenteringRubmaterial material\n";
		return 0;
	}

	return mat;
}

SelfCenteringRub::SelfCenteringRub(int tag, double K1, double K2, double K3,
	double fa, double b, double sd,
	double bd, double rb)
	: UniaxialMaterial(tag, MAT_TAG_SelfCentering),
	k1(K1), k2(K2), k3(K3), ActF(fa), beta(b), SlipDef(sd), BearDef(bd), rBear(rb)
{
	// Find Equivalent Slip Force
	ActDef = ActF / k1;
	SlipF = ActF + (SlipDef - ActDef) * k2;
	if (BearDef != 0) {
		if (SlipDef != 0 && SlipDef < BearDef) {
			BearF = SlipF;
		}
		else {
			BearF = ActF + (BearDef - ActDef) * k2;
		}
	}

	// Initialize variables
	this->revertToStart();
}

SelfCenteringRub::SelfCenteringRub()
	: UniaxialMaterial(0, MAT_TAG_SelfCentering),
	k1(0.0), k2(0.0), k3(0.0), ActF(0.0), beta(0.0),
	SlipDef(0.0), BearDef(0.0), rBear(0.0)
{
	// Initialize variables
	SlipF = 0;
	ActDef = 0;
	BearF = 0;

	this->revertToStart();
}

SelfCenteringRub::~SelfCenteringRub()
{
}

int
SelfCenteringRub::setTrialStrain(double strain, double strainRate)
{
	diffStrain = strain - Cstrain;

	if (fabs(diffStrain) < DBL_EPSILON)
		return 0;

	// Set total strain
	Tstrain = strain;
	noSlipStrain = Tstrain - CslipStrain;

	// Middle Elastic Portion (outside any upper or lower activation)
	//     Entirely elastic response
	if (fabs(noSlipStrain) <= ((1 - beta) * ActF / k1)) {
		Tstress = k1 * noSlipStrain;
		Ttangent = k1;
	}
	else {
		// Positive Quadrant (Top Right) where strain >= 0
		if (noSlipStrain >= 0) {
			// External Fuse Bolt Bearing
			if (BearDef != 0 && Tstrain > BearDef) {
				Tstress = BearF + (Tstrain - BearDef) * rBear * k1;
				Ttangent = rBear * k1;
			}

			// External Fuse Slip
			else if (SlipDef != 0 && noSlipStrain > SlipDef) {
				Tstress = SlipF;
				TslipStrain = CslipStrain + diffStrain;
			}

			// Linear range movement (no upper or
			//     lower activation)
			else if ((noSlipStrain >= ClowerStrainPos) &&
				(noSlipStrain <= CupperStrainPos)) {
				Tstress = (noSlipStrain - CactivStrainPos) * k1;
				Ttangent = k1;
			}
			// Upper Activation
			else if (noSlipStrain > CupperStrainPos) {
				TupperStressPos = ActF +
					(noSlipStrain - ActF/ k1) * k2;
				TupperStrainPos = noSlipStrain;
				TlowerStrainPos = (k1*noSlipStrain-TupperStressPos)/(k1-k3) + (1-beta)*ActF/k1;
				TlowerStressPos = (k1*noSlipStrain-TupperStressPos)*k3/(k1-k3) + (1-beta)*ActF;
				Tstress = TupperStressPos;
				TactivStrainPos = TupperStrainPos - Tstress / k1;

				Ttangent = k2;
			}
			// Lower Activation
			else if (Tstrain < ClowerStrainPos) {
				TlowerStressPos = (1-beta)*ActF +
					(noSlipStrain - (1-beta)*ActF/k1) * k3;
				TlowerStrainPos = noSlipStrain;
				TupperStrainPos = (k1*noSlipStrain-TlowerStressPos)/(k1-k2) + ActF/k1;
				TupperStressPos = (k1*noSlipStrain-TlowerStressPos)*k2/(k1-k2) + ActF;
				Tstress = TlowerStressPos;
				TactivStrainPos = TlowerStrainPos - Tstress / k1;

				Ttangent = k3;
			}
		}

		// Negative Quadrant (Bottom Left) where strain < 0
		else { // Tstrain < 0)
		  // External Fuse Bolt Bearing
			if (BearDef != 0 && Tstrain < -BearDef) {
				Tstress = -BearF + (Tstrain + BearDef) * rBear * k1;
				Ttangent = rBear * k1;
			}

			// External Fuse Slip
			else if (SlipDef != 0 && noSlipStrain < -SlipDef) {
				Tstress = -SlipF;
				TslipStrain = CslipStrain + diffStrain;
			}

			// Linear range movement (no upper or
			//     lower activation)
			else if ((noSlipStrain <= ClowerStrainNeg) &&
				(noSlipStrain >= CupperStrainNeg)) {
				Tstress = (noSlipStrain - CactivStrainNeg) * k1;
				Ttangent = k1;
			}
			// Upper Activation
			else if (noSlipStrain < CupperStrainNeg) {
				TupperStressNeg = -ActF -
					(-noSlipStrain - ActF/k1) * k2;
				TupperStrainNeg = noSlipStrain;
				TlowerStrainNeg = -((k1*(-noSlipStrain)+TupperStressNeg)/(k1-k3) + (1-beta)*ActF/k1);
				TlowerStressNeg = -((k1*(-noSlipStrain)+TupperStressNeg)*k3/(k1-k3) + (1-beta)*ActF);
				Tstress = TupperStressNeg;
				TactivStrainNeg = TupperStrainNeg - Tstress / k1;

				Ttangent = k2;
			}
			// Lower Activation
			else if (Tstrain > ClowerStrainNeg) {
				TlowerStressNeg = -(1-beta)*ActF -
					(-noSlipStrain - (1-beta)*ActF/k1) * k3;
				TlowerStrainNeg = noSlipStrain;
				TupperStrainNeg = -((k1*(-noSlipStrain)+TlowerStressNeg)/(k1-k2) + ActF/k1);
				TupperStressNeg = -((k1*(-noSlipStrain)+TlowerStressNeg)*k2/(k1-k2) + ActF);
				Tstress = TlowerStressNeg;
				TactivStrainNeg = TlowerStrainNeg - Tstress / k1;

				Ttangent = k3;
			}
		}
	}

	return 0;
}

double
SelfCenteringRub::getStress(void)
{
	return Tstress;
}

double
SelfCenteringRub::getTangent(void)
{
	return Ttangent;
}

double
SelfCenteringRub::getStrain(void)
{
	return Tstrain;
}

int
SelfCenteringRub::commitState(void)
{
	// Commit trial history variables
	CactivStrainPos = TactivStrainPos;
	CactivStrainNeg = TactivStrainNeg;
	CslipStrain = TslipStrain;
	CupperStrainPos = TupperStrainPos;
	ClowerStrainPos = TlowerStrainPos;
	CupperStressPos = TupperStressPos;
	ClowerStressPos = TlowerStressPos;
	CupperStrainNeg = TupperStrainNeg;
	ClowerStrainNeg = TlowerStrainNeg;
	CupperStressNeg = TupperStressNeg;
	ClowerStressNeg = TlowerStressNeg;
	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;

	return 0;
}

int
SelfCenteringRub::revertToLastCommit(void)
{
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;

	return 0;
}

int
SelfCenteringRub::revertToStart(void)
{
	// Reset committed history variables
	CactivStrainPos = 0.0;
	CactivStrainNeg = 0.0;
	CslipStrain = 0.0;
	CupperStrainPos = ActDef;
	ClowerStrainPos = (1 - beta) * ActDef;
	CupperStressPos = ActF;
	ClowerStressPos = (1 - beta) * ActF;
	CupperStrainNeg = -CupperStrainPos;
	ClowerStrainNeg = -ClowerStrainPos;
	CupperStressNeg = -CupperStressPos;
	ClowerStressNeg = -ClowerStressPos;

	// Reset trial history variables
	TactivStrainPos = 0.0;
	TactivStrainNeg = 0.0;
	TslipStrain = 0.0;
	TupperStrainPos = ActDef;
	TlowerStrainPos = (1 - beta) * ActDef;
	TupperStressPos = ActF;
	TlowerStressPos = (1 - beta) * ActF;
	TupperStrainNeg = -CupperStrainPos;
	TlowerStrainNeg = -ClowerStrainPos;
	TupperStressNeg = -CupperStressPos;
	TlowerStressNeg = -ClowerStressPos;

	// Initialize state variables
	Tstrain = 0.0;
	Tstress = 0.0;
	Ttangent = k1;

	Cstrain = 0.0;

	return 0;
}

UniaxialMaterial*
SelfCenteringRub::getCopy(void)
{
	SelfCenteringRub* theCopy =
		new SelfCenteringRub(this->getTag(), k1, k2, k3, ActF, beta,
			SlipDef, BearDef, rBear);

	// Copy committed history variables
	theCopy->CactivStrainPos = CactivStrainPos;
	theCopy->CactivStrainNeg = CactivStrainNeg;
	theCopy->CslipStrain = CslipStrain;
	theCopy->CupperStrainPos = CupperStrainPos;
	theCopy->ClowerStrainPos = ClowerStrainPos;
	theCopy->CupperStressPos = CupperStressPos;
	theCopy->ClowerStressPos = ClowerStressPos;
	theCopy->CupperStrainNeg = CupperStrainNeg;
	theCopy->ClowerStrainNeg = ClowerStrainNeg;
	theCopy->CupperStressNeg = CupperStressNeg;
	theCopy->ClowerStressNeg = ClowerStressNeg;

	// Copy trial history variables
	theCopy->TactivStrainPos = TactivStrainPos;
	theCopy->TactivStrainNeg = TactivStrainNeg;
	theCopy->TslipStrain = TslipStrain;
	theCopy->TupperStrainPos = TupperStrainPos;
	theCopy->TlowerStrainPos = TlowerStrainPos;
	theCopy->TupperStressPos = TupperStressPos;
	theCopy->TlowerStressPos = TlowerStressPos;
	theCopy->TupperStrainNeg = TupperStrainNeg;
	theCopy->TlowerStrainNeg = TlowerStrainNeg;
	theCopy->TupperStressNeg = TupperStressNeg;
	theCopy->TlowerStressNeg = TlowerStressNeg;

	// Copy trial state variables
	theCopy->Tstrain = Tstrain;
	theCopy->Tstress = Tstress;
	theCopy->Ttangent = Ttangent;

	theCopy->Cstrain = Cstrain;

	return theCopy;
}

int
SelfCenteringRub::sendSelf(int cTag, Channel& theChannel)
{
	int res = 0;

	static Vector data(27);

	data(0) = this->getTag();
	data(1) = k1;
	data(2) = k2;
	data(3) = k3;
	data(4) = ActF;
	data(5) = beta;
	data(6) = rBear;
	data(7) = SlipDef;
	data(8) = BearDef;
	data(9) = SlipF;
	data(10) = ActDef;
	data(11) = BearF;
	data(12) = CactivStrainPos;
	data(13) = CactivStrainNeg;
	data(14) = CslipStrain;
	data(15) = CupperStrainPos;
	data(16) = ClowerStrainPos;
	data(17) = CupperStressPos;
	data(18) = ClowerStressPos;
	data(19) = CupperStrainNeg;
	data(20) = ClowerStrainNeg;
	data(21) = CupperStressNeg;
	data(22) = ClowerStressNeg;
	data(23) = Tstrain;
	data(24) = Tstress;
	data(25) = Ttangent;
	data(26) = Cstrain;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "SelfCenteringRub::sendSelf() - failed to send data\n";

	return res;
}

int
SelfCenteringRub::recvSelf(int cTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	static Vector data(27);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "SelfCenteringRub::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		this->setTag((int)data(0));
		k1 = data(1);
		k2 = data(2);
		k3 = data(3);
		ActF = data(4);
		beta = data(5);
		rBear = data(6);
		SlipDef = data(7);
		BearDef = data(8);
		SlipF = data(9);
		ActDef = data(10);
		BearF = data(11);
		CactivStrainPos = data(12);
		CactivStrainNeg = data(13);
		CslipStrain = data(14);
		CupperStrainPos = data(15);
		ClowerStrainPos = data(16);
		CupperStressPos = data(17);
		ClowerStressPos = data(18);
		CupperStrainNeg = data(19);
		ClowerStrainNeg = data(20);
		CupperStressNeg = data(21);
		ClowerStressNeg = data(22);
		Tstrain = data(23);
		Tstress = data(24);
		Ttangent = data(25);
		Cstrain = data(26);
	}

	return res;
}

void
SelfCenteringRub::Print(OPS_Stream& s, int flag)
{
	s << "SelfCenteringRub, tag: " << this->getTag() << endln;
	s << "  k1: " << k1 << endln;
	s << "  k2: " << k2 << endln;
	s << "  k3: " << k3 << endln;
	s << "  ActF: " << ActF << endln;
	s << "  beta: " << beta << endln;
	s << "  rBear: " << rBear << endln;
	s << "  SlipDef: " << SlipDef << endln;
	s << "  BearDef: " << BearDef << endln;
}