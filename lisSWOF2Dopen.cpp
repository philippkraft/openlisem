
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-10
#define ve_ca 1e-10

#define GRAV 9.8067
#define EPSILON 1e-10

// 2nd order without iteration dt1, dt2!
double TWorld::fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;

    if (startFlood)
    {
        sumh = getMass(h);

#pragma omp parallel for collapse(2) num_threads(userCores)
        FOR_ROW_COL_MV_L {
            FloodDT->Drc = dt_max;
            FloodT->Drc = 0;
            tma->Drc = 0;
        }

        do {
            // make a copy
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;//std::max(-vmax, std::min(vmax,vx->Drc));
                vys->Drc = vy->Drc;//std::max(-vmax, std::min(vmax,vy->Drc));
            }

            //flow for cells which have h and not done yet (FloodT < _dt)
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (FloodT->Drc < _dt && h->Drc > 0) {
                    double dt = std::max(0.5*FloodDT->Drc, dt_req_min); // use previous timestep to start
                    double vxn, vyn;
                    double vmax = 0.5*(dt < 1 ? _dx+sqrt(_dx/dt) : _dx/dt);  // courant?

                    tma->Drc += 1.0; // nr times a cell is processed

                    //typedef struct vec4 { double v[4]; } vec4;
                    vec4 hll_x1;
                    vec4 hll_x2;
                    vec4 hll_y1;
                    vec4 hll_y2;

                    double dx = ChannelAdj->Drc;
                    double dy = DX->Drc;


                    double H = hs->Drc;
                    double n = N->Drc;
                    double Z = z->Drc;
                    double Vx = std::max(-vmax, std::min(vmax, vxs->Drc));
                    double Vy = std::max(-vmax, std::min(vmax, vys->Drc));

                    bool bc1 = c > 0 && !MV(r,c-1)        ;
                    bool bc2 = c < _nrCols-1 && !MV(r,c+1);
                    bool br1 = r > 0 && !MV(r-1,c)        ;
                    bool br2 = r < _nrRows-1 && !MV(r+1,c);

                    double z_x1 =  bc1 ? z->data[r][c-1] : Z;
                    double z_x2 =  bc2 ? z->data[r][c+1] : Z;
                    double z_y1 =  br1 ? z->data[r-1][c] : Z;
                    double z_y2 =  br2 ? z->data[r+1][c] : Z;

                    double h_x1 =  bc1 ? hs->data[r][c-1] : H;
                    double h_x2 =  bc2 ? hs->data[r][c+1] : H;
                    double h_y1 =  br1 ? hs->data[r-1][c] : H;
                    double h_y2 =  br2 ? hs->data[r+1][c] : H;

                    double vx_x1 = bc1 ? vxs->data[r][c-1] : Vx;
                    double vx_x2 = bc2 ? vxs->data[r][c+1] : Vx;
                    double vx_y1 = br1 ? vxs->data[r-1][c] : Vx;
                    double vx_y2 = br2 ? vxs->data[r+1][c] : Vx;

                    double vy_x1 = bc1 ? vys->data[r][c-1] : Vy;
                    double vy_x2 = bc2 ? vys->data[r][c+1] : Vy;
                    double vy_y1 = br1 ? vys->data[r-1][c] : Vy;
                    double vy_y2 = br2 ? vys->data[r+1][c] : Vy;

                    double fb_x1=0,fb_x2=0,fb_y1=0,fb_y2=0;
                    if (SwitchFlowBarriers) {
                        fb_x1 = std::max(FlowBarrierW->Drc, FlowBarrierE->data[r][c-1]);
                        fb_x2 = std::max(FlowBarrierE->Drc, FlowBarrierE->data[r][c+1]);
                        fb_y1 = std::max(FlowBarrierN->Drc, FlowBarrierS->data[r-1][c]);
                        fb_y2 = std::max(FlowBarrierS->Drc, FlowBarrierN->data[r+1][c]);
                    }

                    vx_x1 = std::max(-vmax, std::min(vmax, vx_x1));
                    vx_x2 = std::max(-vmax, std::min(vmax, vx_x2));
                    vx_y1 = std::max(-vmax, std::min(vmax, vx_y1));
                    vx_y2 = std::max(-vmax, std::min(vmax, vx_y2));

                    vy_x1 = std::max(-vmax, std::min(vmax, vy_x1)); //left
                    vy_x2 = std::max(-vmax, std::min(vmax, vy_x2)); //right
                    vy_y1 = std::max(-vmax, std::min(vmax, vy_y1)); //up
                    vy_y2 = std::max(-vmax, std::min(vmax, vy_y2)); //down

                    // No effect of terrain: use for lakes?
                    // hll_x1 = F_Riemann(h_x1,vx_x1,vy_x1,H,Vx,Vy); // c-1 and c
                    // hll_x2 = F_Riemann(H,Vx,Vy,h_x2,vx_x2,vy_x2); // c and c+1
                    // hll_y1 = F_Riemann(h_y1,vy_y1,vx_y1,H,Vy,Vx); // r-1 and r
                    // hll_y2 = F_Riemann(H,Vy,Vx,h_y2,vy_y2,vx_y2); // r and r+1

                    double fac = 1.0;//DEMdz->Drc; // if Z is in a pit > 10m from the surrounding cells, reduce the effect of the DEM
                    double dz_x1 = fac*(Z - z_x1);
                    double dz_x2 = fac*(z_x2 - Z);
                    double dz_y1 = fac*(Z - z_y1);
                    double dz_y2 = fac*(z_y2 - Z);

                    double h_x1l = std::max(0.0, h_x1 - std::max(0.0,  dz_x1 + fb_x1));
                    double h_x1r = std::max(0.0, H    - std::max(0.0, -dz_x1 + fb_x1));
                    if(bc1)
                        hll_x1 = F_Riemann(h_x1l,vx_x1,vy_x1,h_x1r,Vx,Vy); // c-1 and c
                    else
                        hll_x1 = F_Riemann(0,0,0,h_x1r,Vx,Vy);

                    double h_x2l = std::max(0.0, H    - std::max(0.0,  dz_x2 + fb_x2));
                    double h_x2r = std::max(0.0, h_x2 - std::max(0.0, -dz_x2 + fb_x2));
                    if(bc2)
                        hll_x2 = F_Riemann(h_x2l,Vx,Vy,h_x2r,vx_x2,vy_x2); // c and c+1
                    else
                        hll_x2 = F_Riemann(h_x2l,Vx,Vy,0,0,0);

                    double h_y1u = std::max(0.0, h_y1 - std::max(0.0,  dz_y1 + fb_y1));
                    double h_y1d = std::max(0.0, H    - std::max(0.0, -dz_y1 + fb_y1));
                    if (br1)
                        hll_y1 = F_Riemann(h_y1u,vy_y1,vx_y1,h_y1d,Vy,Vx); // r-1 and r
                    else
                        hll_y1 = F_Riemann(0,0,0,h_y1d,Vy,Vx);

                    double h_y2u = std::max(0.0, H    - std::max(0.0,  dz_y2 + fb_y2));
                    double h_y2d = std::max(0.0, h_y2 - std::max(0.0, -dz_y2 + fb_y2));
                    if(br2)
                        hll_y2 = F_Riemann(h_y2u,Vy,Vx,h_y2d,vy_y2,vx_y2); // r and r+1
                    else
                        hll_y2 = F_Riemann(h_y2u,Vy,Vx,0,0,0);

                    double B = 0.5; //1.0 is theoretical max else faster than gravity
                    double sx_zh_x2 = std::min(B, std::max(-B, (z_x2 + h_x2 - Z - H)/dx));
                    double sy_zh_y1 = std::min(B, std::max(-B, (Z + H - z_y1 - h_y1)/dy));
                    double sx_zh_x1 = std::min(B, std::max(-B, (Z + H - z_x1 - h_x1)/dx));
                    double sy_zh_y2 = std::min(B, std::max(-B, (z_y2 + h_y2 - Z - H)/dy));

                    // if B = 0.5 this can never be >1?
                    double sx_zh = std::min(1.0,std::max(-1.0,limiter(sx_zh_x1, sx_zh_x2)));
                    double sy_zh = std::min(1.0,std::max(-1.0,limiter(sy_zh_y1, sy_zh_y2)));
                    double C = std::min(0.5, courant_factor);
                    double tx = dt/dx;
                    double ty = dt/dy;

                    double flux_x1 = std::max(-H * C,std::min(+tx*hll_x1.v[0],h_x1 * C));
                    double flux_x2 = std::max(-H * C,std::min(-tx*hll_x2.v[0],h_x2 * C));
                    double flux_y1 = std::max(-H * C,std::min(+ty*hll_y1.v[0],h_y1 * C));
                    double flux_y2 = std::max(-H * C,std::min(-ty*hll_y2.v[0],h_y2 * C));

                    double hn = std::max(0.0, H + flux_x1 + flux_x2 + flux_y1 + flux_y2);

                    if(hn > he_ca) {

                        double qxn = H * Vx - tx*(hll_x2.v[1] - hll_x1.v[1]) - ty*(hll_y2.v[2] - hll_y1.v[2])- 0.5 * GRAV *hn*sx_zh * dt;
                        double qyn = H * Vy - tx*(hll_x2.v[2] - hll_x1.v[2]) - ty*(hll_y2.v[1] - hll_y1.v[1])- 0.5 * GRAV *hn*sy_zh * dt;

                        double vsq = sqrt(Vx * Vx + Vy * Vy);
                        double nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01,pow(hn,4.0/3.0));
                        double nsq = nsq1*vsq*dt;

                        vxn = (qxn/(1.0+nsq))/std::max(0.01,hn);
                        vyn = (qyn/(1.0+nsq))/std::max(0.01,hn);


                        if (SwitchTimeavgV) {
                            double fac = 0.5+0.5*std::min(1.0,4*hn)*std::min(1.0,4*hn);
                            fac = fac *exp(- std::max(1.0,dt) / nsq1);
                            vxn = fac * Vx + (1.0-fac) *vxn;
                            vyn = fac * Vy + (1.0-fac) *vyn;
                        }

                        double threshold = 0.01 * _dx; // was 0.01
                        if(hn < threshold)
                        {
                            double h23 = pow(hn, 2.0/3.0);//hn * sqrt(hn)
                            double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
                            double v_kin = (sx_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sx_zh > 0 ? sx_zh : -sx_zh))/(0.001+n);
                            vxn = kinfac * v_kin + vxn*(1.0-kinfac);
                            v_kin = (sy_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sy_zh > 0 ? sy_zh : -sy_zh))/(0.001+n);
                            vyn = kinfac * v_kin + vyn*(1.0-kinfac);
                        }

                    } else {
                        vxn = 0;
                        vyn = 0;
                        hn = 0;
                    }
                    // dan maar even met geweld!
                    if (std::isnan(vxn) || std::isnan(vyn)  )
                    {
                        vxn = 0;
                        vyn = 0;
                        hn= 0;
                    }

                    if (fabs(vxn) <= ve_ca)
                        vxn = 0;
                    if (fabs(vyn) <= ve_ca)
                        vyn = 0;
                    if (hn < he_ca) {
                        vxn = 0;
                        vyn = 0;
                        hn = 0;
                    }

                    // werkt allebij!
                    //double dt_req = courant_factor *_dx/( std::min(dt_max,std::max(0.01,sqrt(vxn*vxn + vyn*vyn))));
                    // gebruik riemann solver cfl
                    double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                    double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]);
                    double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));

                    FloodDT->Drc = dt_req;
                    h->Drc = hn;
                    vx->Drc = vxn;
                    vy->Drc = vyn;
                }
            }


#pragma omp parallel for reduction(min:dt_req_min) collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV {
                double res = FloodDT->Drc;
                dt_req_min = std::min(dt_req_min, res);

            }

            dt_req_min = std::max(TimestepfloodMin, dt_req_min);
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            int cnt = 0;
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (h->Drc > 0) {
                    if (SwitchVariableTimestep)
                        FloodDT->Drc = std::max(dt_req_min,sqrt(dt_req_min*FloodDT->Drc));
                    else
                        FloodDT->Drc = dt_req_min;
                    FloodT->Drc += FloodDT->Drc;
                } else {
                    FloodT->Drc = _dt;
                }
            }

            // nr cells that need processing
#pragma omp parallel for reduction(+:cnt) collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (FloodT->Drc < _dt)
                    cnt++;
            }
            stop = cnt < 1;
                //qDebug() << cnt;

            if (SwitchErosion)
                SWOFSediment(FloodDT,hs,vxs,vys);

            timesum += dt_req_min;
        //     stop = timesum > _dt-0.001;
            count++;

            if(count > F_MaxIter) stop = true;
        } while (!stop);

        correctMassBalance(sumh, h);
    } // if floodstart

    // for screen output
    double avgdt = 0;
#pragma omp parallel for reduction(+:avgdt) collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        FloodDT->Drc = tma->Drc > 0 ? _dt/tma->Drc : dt_max;
        avgdt = avgdt + FloodDT->Drc;
    }
    avgdt = avgdt/nrCells;
    iter_n = count;

    return(avgdt);//_dt/count
}


//              //  if (SwitchMUSCL) {
//                    double dh1   = 0.5*limiter(H-h_x1, h_x2-H);
//                    double dvx1  = 0.5*limiter(Vx-vx_x1, vx_x2-Vx);
//                    double dvy1  = 0.5*limiter(Vy-vy_x1, vy_x2-Vy);

//                    double hlh = H > he_ca ? (H+dh1)/H : 1.0;
//                    vx_x2 = Vx + hlh*dvx1; //xright
//                    vx_x1 = Vx - hlh*dvx1; //xleft
//                    vy_x2 = Vy + hlh*dvy1;
//                    vy_x1 = Vy - hlh*dvy1;

//                    // row -1 and +1
//                    double dh2   = 0.5*limiter(H-h_y1, h_y2-H);
//                    double dvx2  = 0.5*limiter(Vx-vx_y1, vx_y2-Vx);
//                    double dvy2  = 0.5*limiter(Vy-vy_y1, vy_y2-Vy);
//                    double hlh2 = H > he_ca ? (H+dh2)/H : 1.0;
//                    vy_x2 = Vx + hlh2*dvx2;
//                    vy_x1 = Vx - hlh2*dvx2;
//                    vy_y2 = Vy + hlh2*dvy2;
//                    vy_y1 = Vy - hlh2*dvy2;
//              //  }


//                u1r->Drc = _u->Drc + hlh * du;
//                u1l->Drc = _u->Drc - hrh * du;
//                v1r->Drc = _v->Drc + hlh * dv;
//                v1l->Drc = _v->Drc - hrh * dv;
//                h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
//                h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
//                rec = F_Riemann(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);


