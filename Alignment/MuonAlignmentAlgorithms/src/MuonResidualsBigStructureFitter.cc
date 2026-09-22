#ifdef STANDALONE_FITTER
#include "MuonResidualsBigStructureFitter.h"
#else
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsBigStructureFitter.h"
#endif
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "TMath.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

#include <map>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <memory>
#include <fstream>
#include <stdexcept>

#include "Alignment/MuonAlignmentAlgorithms/interface/Tracer.hpp"
#include "Alignment/CommonAlignment/interface/Utilities.h"

#include "TMatrixDSym.h"

static TMinuit* MuonResidualsBigStructureFitter_TMinuit;

static double fcn_zero = 0.0;
static double avg_call_duration = 0.0;
static int call_id = 0;

using DT_6DOF = MuonResidualsBigStructureFitter::DataDT_6DOF;
using DT_5DOF = MuonResidualsBigStructureFitter::DataDT_5DOF;
using CSC_6DOF = MuonResidualsBigStructureFitter::DataCSC_6DOF;
using PARAMS = MuonResidualsBigStructureFitter::PARAMS;

namespace 
{
    TMinuit* minuit;

    // input format: 
    // angle[0] = dphiz; 
    // angle[1] = dphiy;
    // angle[2] = dphix;
    // Reason: rotaions are executed from RIGHT to LEFT (bc of matrix multiplication)
    // first rotate around X by angle[2]
    // than rotate around Y by anlge[1]
    // lastly rotate around Z by angle[0]
    // https://en.wikipedia.org/wiki/Euler_angles#Conversion_to_other_orientation_representations
    // https://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf
    align::RotationType MatrixXYZ(align::EulerAngles const& angles)
    {
        // Scalar is just a typedef for double
        align::Scalar s1 = std::sin(angles[0]), c1 = std::cos(angles[0]);
        align::Scalar s2 = std::sin(angles[1]), c2 = std::cos(angles[1]);
        align::Scalar s3 = std::sin(angles[2]), c3 = std::cos(angles[2]);

        // in constructor(fisrt row of the matrix,     xx xy xz
        //                second row of the matrix,    yx yy yz
        //                third row of the matrix)     zx zy zz
        return align::RotationType(c1*c2, c1*s2*s3 - c3*s1, s1*s3 + c1*c3*s2,
                                c2*s1, c1*c3 + s1*s2*s3, c3*s1*s2 - c1*s3,
                                -s2, c2*s3, c2*c3);
    }

    align::RotationType MatrixXYZ(double dphix, double dphiy, double dphiz)
    {
        align::EulerAngles angles(3);
        angles[0] = dphiz;
        angles[1] = dphiy;
        angles[2] = dphix;

        return MatrixXYZ(angles);
    }

    align::EulerAngles AnglesXYZ(align::RotationType const& matrix)
    {        
        align::EulerAngles res(3);
        double dphix{}; 
        double dphiy{};
        double dphiz{};

        if (matrix.xz() < 1.0)
        {
            if (matrix.xz() > -1.0)
            {
                dphiy = std::asin(-matrix.zx());
                dphiz = std::atan2(matrix.yz(), matrix.xx());
                dphix = std::atan2(matrix.zy(), matrix.zz());
            }
            else 
            {
                Tracer::instance() << "\tWarning in AnglesXYZ at FCN call " << call_id << ": non-unique solution when converting rotation matrix to angles\n";
                Tracer::instance() << matrix << "\n";
                throw std::runtime_error("AnglesXYZ: illegal rotation matrix");
                dphiy = std::asin(1.0)/2;
                dphiz = -std::atan2(-matrix.yz(), matrix.yy());
                dphix = 0.0;
            }
        }
        else 
        {
            Tracer::instance() << "\tWarning in AnglesXYZ at FCN call " << call_id << ": non-unique solution when converting rotation matrix to angles\n";
            Tracer::instance() << matrix << "\n";
            throw std::runtime_error("AnglesXYZ: illegal rotation matrix");
            dphiy = -std::asin(1.0)/2;
            dphiz = std::atan2(-matrix.yz(), matrix.yy());
            dphix = 0.0;
        }

        res[0] = dphix;
        res[1] = dphiy;
        res[2] = dphiz;

        return res;
    }
}

double MuonResidualsBigStructureFitter_logPureGaussian(double residual, double center, double sigma) 
{
    sigma = fabs(sigma);
    static const double cgaus = 0.5 * log(2. * M_PI);
    return (-pow(residual - center, 2) * 0.5 / sigma / sigma) - cgaus - log(sigma);
}

MuonResidualsBigStructureFitter::MuonResidualsBigStructureFitter()
    :   m_DTGeometry(nullptr)
    ,   m_CSCGeometry(nullptr)
    ,   m_printLevel(0)
    ,   m_strategy(0)
    ,   m_value(npar(), 0.0)
    ,   m_error(npar(), 0.0)
    ,   m_loglikelihood(0.0) 
    {}

MuonResidualsBigStructureFitter::MuonResidualsBigStructureFitter(DTGeometry const* DTGeom, 
                                                                 CSCGeometry const* CSCGeom, 
                                                                 std::vector<std::string> const& opt,
                                                                 double min_pt,
                                                                 double max_pt,
                                                                 double errdef)
    :   m_DTGeometry(DTGeom)
    ,   m_CSCGeometry(CSCGeom)
    ,   m_printLevel(0)
    ,   m_strategy(2)
    ,   m_value(npar(), 0.0)
    ,   m_error(npar(), 0.0)
    ,   m_loglikelihood(0.0)
    ,   m_options(opt)
    ,   m_min_pt(min_pt)
    ,   m_max_pt(max_pt)
    ,   m_errdef(errdef)  
    {}

// if something happens comment out m_data, destructor and make it default
MuonResidualsBigStructureFitter::~MuonResidualsBigStructureFitter()
{
    if (!m_data.empty())
    {
        for (auto mapIt: m_data)
        {
            for (auto ptr: mapIt.second)
            {
                delete [] ptr;
            }
        }
    }
}

// NEW FCN START
struct ResidData
{
    double dx;
    double dy;
    double dz;
    double dphix;
    double dphiy;
    double dphiz;
    double track_x;
    double track_y;
    double track_dxdz;
    double track_dydz;
    double slope = 0.0;
    double alpha = 0.0;
};

// functions to compute DT residuals
inline double ResidX(ResidData&& res_data)
{
    auto&& [dx, dy, dz, dphix, dphiy, dphiz, track_x, track_y, track_dxdz, track_dydz, slope, alpha] = res_data;
    return dx - track_dxdz*dz - track_y*track_dxdz*dphix + track_x*track_dxdz*dphiy - track_y*dphiz + slope*alpha;
}

inline double ResidY(ResidData&& res_data)
{
    auto&& [dx, dy, dz, dphix, dphiy, dphiz, track_x, track_y, track_dxdz, track_dydz, slope, alpha] = res_data;
    return dy - track_dydz*dz - track_y*track_dydz*dphix + track_x*track_dydz*dphiy + track_x*dphiz + slope*alpha;
}

inline double ResidDXDZ(ResidData&& res_data)
{
    auto&& [dx, dy, dz, dphix, dphiy, dphiz, track_x, track_y, track_dxdz, track_dydz, slope, alpha] = res_data;
    return -track_dxdz*track_dydz*dphix + (1.0 + track_dxdz*track_dxdz)*dphiy - track_dydz*dphiz;
}

inline double ResidDYDZ(ResidData&& res_data)
{
    auto&& [dx, dy, dz, dphix, dphiy, dphiz, track_x, track_y, track_dxdz, track_dydz, slope, alpha] = res_data;
    return -(1.0 + track_dydz * track_dydz)*dphix + track_dxdz*track_dydz*dphiy + track_dxdz*dphiz;
}

// functions to compute CSC residuals
inline double Resid(ResidData&& res_data, double R)
{
    auto&& [dx, dy, dz, dphix, dphiy, dphiz, track_x, track_y, track_dxdz, track_dydz, slope, alpha] = res_data;
    return dx - (track_x/R - 3.0*(track_x/R)*(track_x/R)*(track_x/R))*dy - track_dxdz*dz - track_y*track_dxdz*dphix + track_x*track_dxdz*dphiy - track_y*dphiz +
            slope*alpha;
}

inline double ResidSlope(ResidData&& res_data, double R)
{
    auto&& [dx, dy, dz, dphix, dphiy, dphiz, track_x, track_y, track_dxdz, track_dydz, slope, alpha] = res_data;
    return -0.5*track_dxdz/R*dy + (track_x/R - track_dxdz*track_dydz)*dphix + (1.0 + track_dxdz*track_dxdz)*dphiy - track_dydz*dphiz;
}

void MuonResidualsBigStructureFitter_FCN(int &npar, double *gin, double &fval, double *par, int iflag) 
{
    std::cout << std::setprecision(8);
    ++call_id;
    MuonResidualsBigStructureFitterFitInfo *fitinfo = (MuonResidualsBigStructureFitterFitInfo *)(minuit->GetObjectFit());
    MuonResidualsBigStructureFitter *gpr_fitter = fitinfo->gpr_fitter();
    DTGeometry const* DTgeom = gpr_fitter->GetDTGeometry();
    CSCGeometry const* CSCgeom = gpr_fitter->GetCSCGeometry();

    fval = 0.0;

    if (!gpr_fitter->HasData())
    {
        return;
    }

    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    for (auto chamber_data = gpr_fitter->DataBegin(); chamber_data != gpr_fitter->DataEnd(); ++chamber_data)
    {
        auto const& [id, resid_vec] = *chamber_data;

        if (id.subdetId() == MuonSubdetId::DT && DTgeom) 
        {
            Surface::PositionType position = DTgeom->idToDet(id)->position();
            Surface::RotationType orientation = DTgeom->idToDet(id)->rotation();

            DTChamberId dtId(id.rawId());
            int station = dtId.station();
            int wheel = dtId.wheel();
                         
            // transform rotation angle from global to local
            // input format: 
            // globAngles[0] = dphiz; 
            // globAngles[1] = dphiy;
            // globAngles[2] = dphix;
            // Reason: rotaions are executed from RIGHT to LEFT (bc of matrix multiplication)
            // first rotate around X by globAngles[2]
            // than rotate around Y by globAngles[1]
            // lastly rotate around Z by globAngles[0]
            
            double dphix = par[static_cast<int>(PARAMS::kAlignPhiX)];
            double dphiy = par[static_cast<int>(PARAMS::kAlignPhiY)];
            double dphiz = par[static_cast<int>(PARAMS::kAlignPhiZ)];

            // align::RotationType glob_rot = MatrixXYZ(globAngles); // proper global gpr matrix in Z(3)Y(2)X(1) convection; 
            align::RotationType glob_rot = MatrixXYZ(dphix, dphiy, dphiz);
            align::RotationType loc_rot = orientation*glob_rot*orientation.transposed();
            align::EulerAngles locAngles = AnglesXYZ(loc_rot);

            // transformed angles
            double alignphix_l = locAngles[0];
            double alignphiy_l = locAngles[1];
            double alignphiz_l = locAngles[2]; 

            // transform translations from global to local
            double dx, dy, dz;
            dx = glob_rot.xx()*position.x() + glob_rot.xy()*position.y() + glob_rot.xz()*position.z() - position.x();
            dy = glob_rot.yx()*position.x() + glob_rot.yy()*position.y() + glob_rot.yz()*position.z() - position.y();
            dz = glob_rot.zx()*position.x() + glob_rot.zy()*position.y() + glob_rot.zz()*position.z() - position.z();

            double alignx_g = par[static_cast<int>(PARAMS::kAlignX)] + dx;
            double aligny_g = par[static_cast<int>(PARAMS::kAlignY)] + dy;
            double alignz_g = par[static_cast<int>(PARAMS::kAlignZ)] + dz;

            GlobalVector translation_g = GlobalVector(alignx_g, aligny_g, alignz_g);
            LocalVector translation_l = DTgeom->idToDet(id)->toLocal(translation_g);

            // transformed shifts
            double alignx_l = translation_l.x(); 
            double aligny_l = translation_l.y(); 
            double alignz_l = translation_l.z();

            if (station != 4)
            {
                std::vector<double> sigmas = gpr_fitter->GetStdDevs(id);
                double resXsigma = sigmas[static_cast<int>(DT_6DOF::kResidX)];
                double resYsigma = sigmas[static_cast<int>(DT_6DOF::kResidY)];
                double slopeXsigma = sigmas[static_cast<int>(DT_6DOF::kResSlopeX)];
                double slopeYsigma = sigmas[static_cast<int>(DT_6DOF::kResSlopeY)];

                std::vector<double*>::const_iterator resid_begin = resid_vec.cbegin();
                std::vector<double*>::const_iterator resid_end = resid_vec.cend();
                std::vector<double*>::const_iterator it = resid_begin;
                for (it = resid_begin; it != resid_end; ++it)
                {
                    double residX = (*it)[static_cast<int>(DT_6DOF::kResidX)];
                    double residY = (*it)[static_cast<int>(DT_6DOF::kResidY)];
                    double resslopeX = (*it)[static_cast<int>(DT_6DOF::kResSlopeX)];
                    double resslopeY = (*it)[static_cast<int>(DT_6DOF::kResSlopeY)];
                    double positionX = (*it)[static_cast<int>(DT_6DOF::kPositionX)];
                    double positionY = (*it)[static_cast<int>(DT_6DOF::kPositionY)];
                    double angleX = (*it)[static_cast<int>(DT_6DOF::kAngleX)];
                    double angleY = (*it)[static_cast<int>(DT_6DOF::kAngleY)];

                    double residXpeak = ResidX({alignx_l, aligny_l, alignz_l, alignphix_l, alignphiy_l, alignphiz_l, positionX, positionY, angleX, angleY, resslopeX});
                    double residYpeak = ResidY({alignx_l, aligny_l, alignz_l, alignphix_l, alignphiy_l, alignphiz_l, positionX, positionY, angleX, angleY, resslopeY});
                    double slopeXpeak = ResidDXDZ({alignx_l, aligny_l, alignz_l, alignphix_l, alignphiy_l, alignphiz_l, positionX, positionY, angleX, angleY});
                    double slopeYpeak = ResidDYDZ({alignx_l, aligny_l, alignz_l, alignphix_l, alignphiy_l, alignphiz_l, positionX, positionY, angleX, angleY});

                    double weight = 1.0;

                    fval += -weight * MuonResidualsBigStructureFitter_logPureGaussian(residX, residXpeak, resXsigma);
                    fval += -weight * MuonResidualsBigStructureFitter_logPureGaussian(residY, residYpeak, resYsigma);
                    fval += -weight * MuonResidualsBigStructureFitter_logPureGaussian(resslopeX, slopeXpeak, slopeXsigma);
                    fval += -weight * MuonResidualsBigStructureFitter_logPureGaussian(resslopeY, slopeYpeak, slopeYsigma);
                }
            }
            else 
            {
                continue;
            }
        }

        if (id.subdetId() == MuonSubdetId::CSC && CSCgeom)
        {
            Surface::PositionType position = CSCgeom->idToDet(id)->position();
            Surface::RotationType orientation = CSCgeom->idToDet(id)->rotation();

            double csc_x = position.x();
            double csc_y = position.y();

            double csc_R = sqrt(csc_x*csc_x + csc_y*csc_y);

            // transform rotation angle from global to local
            double dphix = par[static_cast<int>(PARAMS::kAlignPhiX)];
            double dphiy = par[static_cast<int>(PARAMS::kAlignPhiY)];
            double dphiz = par[static_cast<int>(PARAMS::kAlignPhiZ)];

            // align::RotationType glob_rot = MatrixXYZ(globAngles); // proper global gpr matrix in Z(3)Y(2)X(1) convection; 
            align::RotationType glob_rot = MatrixXYZ(dphix, dphiy, dphiz);
            align::RotationType loc_rot = orientation*glob_rot*orientation.transposed();
            align::EulerAngles locAngles = AnglesXYZ(loc_rot);

            // transformed angles
            double alignphix_l = locAngles[0];
            double alignphiy_l = locAngles[1];
            double alignphiz_l = locAngles[2]; 

            // transform translations from global to local
            double dx, dy, dz;
            dx = glob_rot.xx()*position.x() + glob_rot.xy()*position.y() + glob_rot.xz()*position.z() - position.x();
            dy = glob_rot.yx()*position.x() + glob_rot.yy()*position.y() + glob_rot.yz()*position.z() - position.y();
            dz = glob_rot.zx()*position.x() + glob_rot.zy()*position.y() + glob_rot.zz()*position.z() - position.z();

            double alignx_g = par[static_cast<int>(PARAMS::kAlignX)] + dx;
            double aligny_g = par[static_cast<int>(PARAMS::kAlignY)] + dy;
            double alignz_g = par[static_cast<int>(PARAMS::kAlignZ)] + dz;

            GlobalVector translation_g = GlobalVector(alignx_g, aligny_g, alignz_g);
            LocalVector translation_l = CSCgeom->idToDet(id)->toLocal(translation_g);

            // transformed shifts
            double alignx_l = translation_l.x(); 
            double aligny_l = translation_l.y(); 
            double alignz_l = translation_l.z();

            std::vector<double> sigmas = gpr_fitter->GetStdDevs(id);
            double resSigma = sigmas[static_cast<int>(CSC_6DOF::kResid)];
            double resSlopeSigma = sigmas[static_cast<int>(CSC_6DOF::kResSlope)];
            
            std::vector<double*>::const_iterator resid_begin = resid_vec.cbegin();
            std::vector<double*>::const_iterator resid_end = resid_vec.cend();
            std::vector<double*>::const_iterator it = resid_begin;
            for (it = resid_begin; it != resid_end; ++it)
            {
                double resid = (*it)[static_cast<int>(CSC_6DOF::kResid)];
                double resSlope = (*it)[static_cast<int>(CSC_6DOF::kResSlope)];
                double positionX = (*it)[static_cast<int>(CSC_6DOF::kPositionX)];
                double positionY = (*it)[static_cast<int>(CSC_6DOF::kPositionY)];
                double angleX = (*it)[static_cast<int>(CSC_6DOF::kAngleX)];
                double angleY = (*it)[static_cast<int>(CSC_6DOF::kAngleY)];

                double residPeak = Resid({alignx_l, aligny_l, alignz_l, alignphix_l, alignphiy_l, alignphiz_l, positionX, positionY, angleX, angleY}, csc_R);
                double resSlopePeak = ResidSlope({alignx_l, aligny_l, alignz_l, alignphix_l, alignphiy_l, alignphiz_l, positionX, positionY, angleX, angleY}, csc_R);

                double weight = 1.0;

                fval += -weight * MuonResidualsBigStructureFitter_logPureGaussian(resid, residPeak, resSigma);
                fval += -weight * MuonResidualsBigStructureFitter_logPureGaussian(resSlope, resSlopePeak, resSlopeSigma);
            }
        }  
    }
    
    if (call_id > 1) fval -= fcn_zero;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    std::cout << "******************************************" << "\n";
    std::cout << "FCN call #: " << call_id << "\n";
    std::cout << "execution time " << duration.count() << " ms" << "\n";
    std::cout << "FCN = " << fval << "\n";
    std::cout << "global params: ";
    for (int i = 0; i < npar; ++i)
    {
        std::cout << par[i] << " ";
    } 
    std::cout << "\n";
    std::cout << "******************************************" << "\n";
    avg_call_duration += static_cast<double>(duration.count());
}
// NEW FCN END

void MuonResidualsBigStructureFitter::inform(TMinuit* tMinuit) { minuit = tMinuit; }

bool MuonResidualsBigStructureFitter::dofit(void (*fcn)(int &, double *, double &, double *, int),
                                   std::vector<int> &parNum,
                                   std::vector<std::string> &parName,
                                   std::vector<double> &start,
                                   std::vector<double> &step,
                                   std::vector<double> &low,
                                   std::vector<double> &high) 
{
    MuonResidualsBigStructureFitterFitInfo *fitinfo = new MuonResidualsBigStructureFitterFitInfo(this);

    // configure Minuit object
    MuonResidualsBigStructureFitter_TMinuit = new TMinuit(npar());
    MuonResidualsBigStructureFitter_TMinuit->SetPrintLevel();
    MuonResidualsBigStructureFitter_TMinuit->SetObjectFit(fitinfo);
    MuonResidualsBigStructureFitter_TMinuit->SetFCN(fcn);
    inform(MuonResidualsBigStructureFitter_TMinuit);

    std::vector<int>::const_iterator iNum = parNum.begin();
    std::vector<std::string>::const_iterator iName = parName.begin();
    std::vector<double>::const_iterator istart = start.begin();
    std::vector<double>::const_iterator istep = step.begin();
    std::vector<double>::const_iterator ilow = low.begin();
    std::vector<double>::const_iterator ihigh = high.begin();

    for (; iNum != parNum.end(); ++iNum, ++iName, ++istart, ++istep, ++ilow, ++ihigh) 
    {
        MuonResidualsBigStructureFitter_TMinuit->DefineParameter(*iNum, iName->c_str(), *istart, *istep, *ilow, *ihigh);
        if (*iNum > 5)
        {
            MuonResidualsBigStructureFitter_TMinuit->FixParameter(*iNum);
        }
    }

    Tracer::instance() << "start: ";
    Tracer::instance().copy(start.begin(), start.end());
    Tracer::instance() << "\n";

    Double_t fmin, fedm, errdef;
    Int_t npari, nparx, istat;

    double arglist[10];
    int ierflg;
    int smierflg;  //second MIGRAD ierflg

    // chi^2 errors should be 1.0, log-likelihood should be 0.5
    for (int i = 0; i < 10; i++)
    {
        arglist[i] = 0.0;
    }
    arglist[0] = m_errdef; // 68% CL with one parameter <=> 0.5 
    ierflg = 0;
    smierflg = 0;
    MuonResidualsBigStructureFitter_TMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
    if (ierflg != 0) 
    {
        delete MuonResidualsBigStructureFitter_TMinuit;
        delete fitinfo;
        return false;
    }

    // set strategy = 2 (more refined fits)
    for (int i = 0; i < 10; i++)
    {
        arglist[i] = 0.0;
    }
    arglist[0] = m_strategy;
    ierflg = 0;
    MuonResidualsBigStructureFitter_TMinuit->mnexcm("SET STR", arglist, 1, ierflg);
    if (ierflg != 0) 
    {
        delete MuonResidualsBigStructureFitter_TMinuit;
        delete fitinfo;
        return false;
    }

    int ndim = npar();
    double par[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int flag = 1;
    double grad[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MuonResidualsBigStructureFitter_TMinuit->Eval(ndim, grad, fcn_zero, par, flag);
    Tracer::instance() << std::setprecision(10) << "FCN(0) = " << fcn_zero << "\n";

    bool try_again = false;

    // minimize
    for (int i = 0; i < 10; i++)
    {
        arglist[i] = 0.0;
    }
    arglist[0] = 50000;
    ierflg = 0;
    MuonResidualsBigStructureFitter_TMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
    MuonResidualsBigStructureFitter_TMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);

    Tracer::instance() << "MIGRAD:\n"
                       << "\tfmin = " << fmin << "\n"
                       << "\tfedm = " << fedm << "\n"
                       << "\tistat = " << istat << "\n"
                       << "\tierflg = " << ierflg << "\n";
                       
    if (ierflg != 0)
    {
        try_again = true;
    }
    // just once more, if needed (using the final Minuit parameters from the failed fit; often works)
    if (try_again) 
    {
        for (int i = 0; i < 10; i++)
        {
            arglist[i] = 0.0;
        }
        arglist[0] = 50000;
        arglist[1] = 0.05;
        smierflg = 0;
        MuonResidualsBigStructureFitter_TMinuit->mnexcm("MIGRAD", arglist, 2, smierflg);
    }

    for (int i = 0; i < 10; i++)
    {
        arglist[i] = 0.0;
    }
    ierflg = 0;
    MuonResidualsBigStructureFitter_TMinuit->mnexcm("HESSE", arglist, 0, ierflg);
    MuonResidualsBigStructureFitter_TMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);

    Tracer::instance() << "HESSE:\n"
                       << "\tfmin = " << fmin << "\n"
                       << "\tfedm = " << fedm << "\n"
                       << "\tistat = " << istat << "\n"
                       << "\tierflg = " << ierflg << "\n";

    // read-out the results
    m_loglikelihood = -fmin;

    m_value.clear();
    m_error.clear();
    for (int i = 0; i < npar(); i++) 
    {
        double v, e;
        MuonResidualsBigStructureFitter_TMinuit->GetParameter(i, v, e);
        m_value.push_back(v);
        m_error.push_back(e);
    }

    Tracer::instance() << "m_value: ";
    Tracer::instance().copy(m_value.begin(), m_value.end());
    Tracer::instance() << "\n";
    Tracer::instance() << "m_error: ";
    Tracer::instance().copy(m_error.begin(), m_error.end());
    Tracer::instance() << "\n";

    avg_call_duration /= call_id;
    Tracer::instance() << "average call duration = " << avg_call_duration << " ms\n";
    Tracer::instance() << "# FCN calls = " << call_id << "\n";

    delete MuonResidualsBigStructureFitter_TMinuit;
    delete fitinfo;
    if (smierflg != 0)
    {
        return false;
    }    
    return true;
}

bool MuonResidualsBigStructureFitter::Fit()
{
    std::vector<int> nums{ static_cast<int>(PARAMS::kAlignX),
                           static_cast<int>(PARAMS::kAlignY),
                           static_cast<int>(PARAMS::kAlignZ),
                           static_cast<int>(PARAMS::kAlignPhiX),
                           static_cast<int>(PARAMS::kAlignPhiY),
                           static_cast<int>(PARAMS::kAlignPhiZ) };

    std::vector<std::string> names{ "AlignX",
                                    "AlignY",
                                    "AlignZ",
                                    "AlignPhiX",
                                    "AlignPhiY",
                                    "AlignPhiZ" };

    std::vector<double> starts{ 0.01,
                                0.01,
                                0.01,
                                0.0001,
                                0.0001,
                                0.0001 };

    std::vector<double> steps{ 0.002,
                               0.002,
                               0.002,
                               0.00002,
                               0.00002,
                               0.00002 };

    // parameter ranges 
    std::vector<double> lows{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::vector<double> highs{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    return dofit(&MuonResidualsBigStructureFitter_FCN, nums, names, starts, steps, lows, highs);
}

size_t MuonResidualsBigStructureFitter::NTypesOfResid(Alignable const* ali) const
{
    DetId id = ali->geomDetId();
    if (id.subdetId() == MuonSubdetId::CSC) return static_cast<size_t>(DataCSC_6DOF::kNData);
    if (id.subdetId() == MuonSubdetId::DT)
    {
        DTChamberId dtId(id.rawId());
        if (dtId.station() != 4) return static_cast<size_t>(DataDT_6DOF::kNData);
        return static_cast<size_t>(DataDT_5DOF::kNData);
    }
    return 0;
}

size_t MuonResidualsBigStructureFitter::NTypesOfResid(DetId const& id) const
{
    if (id.subdetId() == MuonSubdetId::CSC) return static_cast<size_t>(DataCSC_6DOF::kNData);
    if (id.subdetId() == MuonSubdetId::DT)
    {
        DTChamberId dtId(id.rawId());
        if (dtId.station() != 4) return static_cast<size_t>(DataDT_6DOF::kNData);
        return static_cast<size_t>(DataDT_5DOF::kNData);
    }
    return 0;
}

void MuonResidualsBigStructureFitter::CopyData(std::map<Alignable*, MuonResidualsTwoBin*> const& from)
{
    for (auto fitIt = from.cbegin(); fitIt != from.cend(); ++fitIt)
    {
        // std::vector<double*> contains many arrays, each of size NTypesOfResid(ali); each array holds residuals calculated from one track
        // in other words, size of each double* is number of different residuals in a chamber (depends on chamber type)
        // size of std::vector<double*> - number of tracks in the chamber from which residuals were calculated

        Alignable* ali = fitIt->first;
        DetId id = ali->geomDetId();

        // don't copy what's not needed
        if (!IsSelected(id))
        {
            continue;
        }

        // sizes of arrays in vector
        size_t resPosSz = fitIt->second->numResidualsPos();
        size_t resNegSz = fitIt->second->numResidualsNeg();

        // select only fraction of residuals
        size_t totResSz = resPosSz + resNegSz;

        if (totResSz <= 2)
        {
            continue;
        }

        int pt_idx = -1;
        if (id.subdetId() == MuonSubdetId::CSC)
        {
            pt_idx = MuonResiduals5DOFFitter::kPt;
        }
        else if (id.subdetId() == MuonSubdetId::DT)
        {
            DTChamberId dtId(id.rawId());
            if (dtId.station() != 4)
            {
                pt_idx = MuonResiduals6DOFFitter::kPt;
            }
            else
            {
                pt_idx = MuonResiduals5DOFFitter::kPt;
            }
        }
        
        if (pt_idx < 0)
        {
            throw std::runtime_error("MuonResidualsBigStructureFitter::CopyData: invalid index of pt in array of alignment data.");
        }

        // number type
        size_t numResTypes = NTypesOfResid(ali);
        std::vector<double*> tmp;
        tmp.reserve(totResSz);
        auto resPosBegin = fitIt->second->residualsPos_begin();
        auto resNegBegin = fitIt->second->residualsNeg_begin();
        for (size_t idx = 0; idx < totResSz; ++idx)
        {
            double* arr = *(resPosBegin + idx);
            double pt = arr[pt_idx];

            if (pt < m_min_pt || pt > m_max_pt)
            {
                continue;
            }

            double* data = new double[numResTypes];
            if (idx < resPosSz) 
            {
                std::copy(*(resPosBegin + idx), *(resPosBegin + idx) + numResTypes, data);
            }
            else
            {
                size_t negResIdx = idx - resPosSz;
                std::copy(*(resNegBegin + negResIdx), *(resNegBegin + negResIdx) + numResTypes, data);
            }
            tmp.push_back(data);
        } 
        m_data.insert({id, tmp});
    }
}

int MuonResidualsBigStructureFitter::TrackCount() const
{
    int cnt = 0;
    for (auto const& [id, resid_vec]: m_data)
    {
        cnt += resid_vec.size();
    }
    return cnt;
}

void MuonResidualsBigStructureFitter::ReleaseData()
{
    for (auto& [id, resid_vec]: m_data)
    {
        for (auto res_ptr: resid_vec)
        {
            delete [] res_ptr;
        }
    }

    m_data.clear();
}

bool MuonResidualsBigStructureFitter::IsSelected(DetId const& id) const
{
    if (m_options.empty()) return false;

    if (id.subdetId() == MuonSubdetId::CSC)
    {
        bool ringSelec, statSelec, endcSelec;
        CSCDetId cscId(id.rawId());
        int s = cscId.station();
        int e = cscId.endcap();
        int r = cscId.ring();
        endcSelec = m_options[Options::CSCEndcaps].empty() ? true : (m_options[Options::CSCEndcaps][e - 1] == '1');
        ringSelec = m_options[Options::CSCRings].empty() ? true : (m_options[Options::CSCRings][r - 1] == '1');
        statSelec = m_options[Options::CSCStations].empty() ? true : (m_options[Options::CSCStations][s - 1] == '1');
        return (endcSelec && ringSelec && statSelec); 
    }
    if (id.subdetId() == MuonSubdetId::DT)
    {
        DTChamberId dtId(id.rawId());
        int s = dtId.station();
        int w = dtId.wheel();
        // -2 <= w <= 2
        bool wheelSelec, statSelec;
        // if wheel string is empty, allow all wheels; same for stations
        wheelSelec = m_options[Options::DTWheels].empty() ? true : (m_options[Options::DTWheels][w + 2] == '1');
        statSelec = m_options[Options::DTStations].empty() ? true : (m_options[Options::DTStations][s - 1] == '1');
        return (wheelSelec && statSelec); 
    }
    return false;
}

void MuonResidualsBigStructureFitter::SetOptions(std::vector<std::string> const& options)
{
    if (options.size() != static_cast<size_t>(Options::OptCount))
    {
        throw std::runtime_error("Wrong number of options provided");
    }

    m_options = options;
}

void MuonResidualsBigStructureFitter::SetOption(size_t optId, std::string const& option)
{
    if (optId >= static_cast<size_t>(Options::OptCount))
    {
        throw std::runtime_error("Attempted to set non-existent option");
    } 
    if (m_options.empty())
    {
        m_options.resize(static_cast<size_t>(Options::OptCount));
    }
    m_options[optId] = option;
}

void MuonResidualsBigStructureFitter::CalcStats()
{
    if (m_data.empty()) return;

    double val{};

    for (auto const& [id, resid_vec]: m_data)
    {
        if (resid_vec.empty())
        {
            continue;
        }

        int n = resid_vec.size();

        if (id.subdetId() == MuonSubdetId::DT)
        {
            DTChamberId dtId(id.rawId());
            int station = dtId.station();

            if (station != 4)
            {
                double sum_sqr_dx = 0.0;
                double sum_dx = 0.0;

                double sum_sqr_dy = 0.0;
                double sum_dy = 0.0;

                double sum_sqr_dxdz = 0.0;
                double sum_dxdz = 0.0;

                double sum_sqr_dydz = 0.0;
                double sum_dydz = 0.0;

                for (auto res_ptr: resid_vec)
                {
                    val = res_ptr[static_cast<int>(DT_6DOF::kResidX)];
                    sum_sqr_dx += val*val;
                    sum_dx += val;

                    val = res_ptr[static_cast<int>(DT_6DOF::kResidY)];
                    sum_sqr_dy += val*val;
                    sum_dy += val;

                    val = res_ptr[static_cast<int>(DT_6DOF::kResSlopeX)];
                    sum_sqr_dxdz += val*val;
                    sum_dxdz += val;

                    val = res_ptr[static_cast<int>(DT_6DOF::kResSlopeY)];
                    sum_sqr_dydz += val*val;
                    sum_dydz += val;
                }

                double residX_mean = sum_dx/n;
                double residY_mean = sum_dy/n;
                double resslopeX_mean = sum_dxdz/n;
                double resslopeY_mean = sum_dydz/n;

                double residX_stddev = sqrt(sum_sqr_dx/n - sum_dx*sum_dx/(n*n));
                double residY_stddev = sqrt(sum_sqr_dy/n - sum_dy*sum_dy/(n*n));
                double resslopeX_stddev = sqrt(sum_sqr_dxdz/n - sum_dxdz*sum_dxdz/(n*n));
                double resslopeY_stddev = sqrt(sum_sqr_dydz/n - sum_dydz*sum_dydz/(n*n));

                std::vector<double> stds{residX_stddev, residY_stddev, resslopeX_stddev, resslopeY_stddev};
                std::vector<double> means{residX_mean, residY_mean, resslopeX_mean, resslopeY_mean};
                m_stats.insert({id, {std::move(means), std::move(stds)}});
            }
            else
            {
                double sum_sqr_dx = 0.0;
                double sum_dx = 0.0;

                double sum_sqr_dxdz = 0.0;
                double sum_dxdz = 0.0;

                for (auto res_ptr: resid_vec)
                {
                    val = res_ptr[static_cast<int>(DT_5DOF::kResid)];
                    sum_sqr_dx += val*val;
                    sum_dx += val;

                    val = res_ptr[static_cast<int>(DT_5DOF::kResSlope)];
                    sum_sqr_dxdz += val*val;
                    sum_dxdz += val;
                }

                double resid_mean = sum_dx/n;
                double resslope_mean = sum_dxdz/n;

                double resid_stddev = sqrt(sum_sqr_dx/n - sum_dx*sum_dx/(n*n));
                double resslope_stddev = sqrt(sum_sqr_dxdz/n - sum_dxdz*sum_dxdz/(n*n));

                std::vector<double> stds{resid_stddev, resslope_stddev};
                std::vector<double> means{resid_mean, resslope_mean};
                m_stats.insert({id, {std::move(means), std::move(stds)}});
            }
        }

        if (id.subdetId() == MuonSubdetId::CSC)
        {
            double sum_sqr_dx = 0.0;
            double sum_dx = 0.0;

            double sum_sqr_dxdz = 0.0;
            double sum_dxdz = 0.0;

            for (auto res_ptr: resid_vec)
            {
                val = res_ptr[static_cast<int>(CSC_6DOF::kResid)];
                sum_sqr_dx += val*val;
                sum_dx += val;

                val = res_ptr[static_cast<int>(CSC_6DOF::kResSlope)];
                sum_sqr_dxdz += val*val;
                sum_dxdz += val;
            }

            double resid_mean = sum_dx/n;
            double resslope_mean = sum_dxdz/n;

            double resid_stddev = sqrt(sum_sqr_dx/n - sum_dx*sum_dx/(n*n));
            double resslope_stddev = sqrt(sum_sqr_dxdz/n - sum_dxdz*sum_dxdz/(n*n));

            std::vector<double> stds{resid_stddev, resslope_stddev};
            std::vector<double> means{resid_mean, resslope_mean};
            m_stats.insert({id, {std::move(means), std::move(stds)}});
        }
    }
}
