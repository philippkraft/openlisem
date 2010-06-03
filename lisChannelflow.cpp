/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
Functionality in lisChannelflow.cpp:
- fraction of water and sediment flowing into the channel
- calc V, alpha and Q in the channel
- calc channelflow, channelheight, kin wave
*/

#include "model.h"

//---------------------------------------------------------------------------
//fraction of water and sediment flowing into the channel
void TWorld::ToChannel(void)
{
	if (SwitchIncludeChannel)
	{
		FOR_ROW_COL_MV_CH
		{
			double fractiontochannel = min(_dt*V->Drc/(0.5*(_dx-ChannelWidthUpDX->Drc)), 1.0);
			double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

			if (SwitchAllinChannel)
				if (Outlet->Drc == 1)
					fractiontochannel = 1.0;
			// in catchment outlet cell, throw everything in channel

			if (SwitchBuffers)
				if (BufferID->Drc > 0)
					fractiontochannel = 1.0;
			// where there is a buffer in the channel, all goes in the channel

			RunoffVolinToChannel->Drc = fractiontochannel*Volume;
			// water diverted to the channel
			WHrunoff->Drc *= (1-fractiontochannel);
			// adjust water height
			if (SwitchErosion)
			{
				SedToChannel->Drc = fractiontochannel*Sed->Drc;
				//sediment diverted to the channel
				Sed->Drc -= SedToChannel->Drc;
				// adjust sediment in suspension
			}
		}
		CalcVelDisch();
		// recalc velocity and discharge
	}
}
//---------------------------------------------------------------------------
// V, alpha and Q in the channel
void TWorld::CalcVelDischChannel()
{
	/*
    dw      FW      dw
   \  |            |  /
    \ |         wh | /
     \|____________|/
	 */
	FOR_ROW_COL_MV_CH
	{
		double Perim, Radius, Area, beta = 0.6;
		double _23 = 2.0/3.0;
		double beta1 = 1/beta;
		double wh = ChannelWH->Drc;
		double FW = ChannelWidth->Drc;
		double grad = sqrt(ChannelGrad->Drc);
		double dw = 0.5*(ChannelWidthUpDX->Drc - FW); // extra width when non-rectamgular

		if (dw > 0)
		{
			Perim = FW + 2*sqrt(wh*wh + dw*dw);
			Area = FW*wh + wh*dw*2;
		}
		else
		{
			Perim = FW + 2*wh;
			Area = FW*wh;
		}

		//Perim = ChannelWidth->Drc + 2*wh/cos(atan(ChannelSide->Drc));
		// cos atanb more expensive than sqrt ?
		//Area = ChannelWidth->Drc*wh + wh*(ChannelWidthUpDX->Drc - ChannelWidth->Drc);


		if (Perim > 0)
			Radius = Area/Perim;
		else
			Radius = 0;

		ChannelAlpha->Drc = pow(ChannelN->Drc/grad * powl(Perim, _23),beta);

		if (ChannelAlpha->Drc > 0)
			ChannelQ->Drc = pow(Area/ChannelAlpha->Drc, beta1);
		else
			ChannelQ->Drc = 0;

		ChannelV->Drc = pow(Radius, _23)*grad/ChannelN->Drc;
	}
	else
	{
		ChannelAlpha->Drc = 0;
		ChannelQ->Drc = 0;
		ChannelV->Drc = 0;
	}
}
//---------------------------------------------------------------------------
//- calc channelflow, channelheight, kin wave
void TWorld::ChannelFlow(void)
{
	if (!SwitchIncludeChannel)
		return;

	// calculate new channel WH , WidthUp and Volume
	FOR_ROW_COL_MV_CH
	{
		/*---- Water ----*/

		ChannelQsn->Drc =0;
		Channelq->Drc =0;
		ChannelQoutflow->Drc =0;
		ChannelWH->Drc = 0;

		ChannelWaterVol->Data[r][c] += RunoffVolinToChannel->Drc;
		// add inflow to channel
		if (ChannelWidth->Drc > 0)
			ChannelWaterVol->Drc += Rainc->Drc*ChannelWidthUpDX->Drc*DX->Drc;
		// add rainfall in m3, no interception

		if (ChannelSide->Drc == 0 && ChannelWidth->Drc > 0)// rectangular channel
		{
			ChannelWidthUpDX->Drc = ChannelWidth->Drc;
			ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*DX->Drc);
		}
		else  // non-rectangular
		{
			if (ChannelWaterVol->Drc > 0)
			{
				double a = ChannelSide->Drc*DX->Drc/ChannelWaterVol->Drc;
				double b = ChannelWidth->Drc*DX->Drc/ChannelWaterVol->Drc;
				double cc = -1.0;
				if (a > 0)
					ChannelWH->Drc = -b + sqrt(b*b-4*a*-cc)/(2*a);
				else
					ChannelWH->Drc = 0;
			}
			// new WH with abc method
		}

		if (ChannelWidth->Drc > 0)
			ChannelWidthUpDX->Drc = min(0.9*_dx, ChannelWidth->Drc+2*ChannelSide->Drc*ChannelWH->Drc);
		// new channel width with new WH, goniometric, side is top angle tan, 1 is 45 degr
		// cannot be more than 0.9*_dx

		if (RoadWidthDX->Drc > 0)
			ChannelWidthUpDX->Drc = min(0.9*_dx-RoadWidthDX->Drc, ChannelWidthUpDX->Drc);
		// channel cannot be wider than _dx-road
		//TODO zit al in gridcell, nodig hier?

	}

	CalcVelDischChannel();

	/*---- Sediment ----*/

	if (SwitchErosion)
	{
		ChannelFlowDetachment();

		FOR_ROW_COL_MV_CH
		{
			ChannelQsoutflow->Drc = 0;
			ChannelQs->Drc = ChannelQ->Drc * ChannelConc->Drc;
		}
	}

	ChannelQn->setMV();

	FOR_ROW_COL_MV_CH
	{
		if (LDDChannel->Drc == 5)
		{
			Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, DX,
					ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed);

			ChannelQoutflow->Drc = ChannelQn->Drc * _dt;
        	if (SwitchErosion)
        		ChannelQsoutflow->Drc = ChannelQsn->Drc * _dt;
			// these maps now contain m3 and kg per timestep in pit cells
		}
	}
	ChannelQn->cover(0); // avoid missing values around channel for adding to Qn for output
	ChannelQs->cover(0);

	FOR_ROW_COL_MV_CH
	{
		double ChannelArea = ChannelAlpha->Drc*pow(ChannelQn->Drc, 0.6);
		ChannelWH->Drc = ChannelArea/ChannelWidthUpDX->Drc;
		// water height is not used except for output e.g. watervolume is cycled

		ChannelWaterVol->Drc = ChannelArea * DX->Drc;
		// total water vol after kin wave in m3, going to the next timestep

		if (SwitchErosion)
		{
			ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc, ChannelDep->Drc);
			// correct for very high concentrations, max 850 g/l
		}
	}
}
//---------------------------------------------------------------------------
