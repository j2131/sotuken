#pragma once

#include <vector>

namespace {
	constexpr PS::F64 kBoxHalfX = 0.002;  // 0.2 cm
	constexpr PS::F64 kBoxHalfY = 0.025;  // 2.5 cm (length 5 cm)
	//constexpr PS::F64 kBoxHalfZ = 0.01;  // 1.0 cm
    constexpr PS::F64 kBoxHalfZ = 0.02;  // 2.0 cm
	constexpr PS::F64 kVibrationOmega = 2.0 * M_PI;
	//constexpr PS::F64vec kVibrationAmp(0.001, 0.0005, 0.001);
    //constexpr PS::F64vec kVibrationAmp(0.0, 0.0, kBoxHalfZ*0.1);
    constexpr PS::F64vec kVibrationAmp(0.0, 0.0, 0.0);
    constexpr PS::F64 kBoxVolume = (2.0 * kBoxHalfX) * (2.0 * kBoxHalfY) * (2.0 * kBoxHalfZ);
}

template <class ThisPtcl> class Problem{
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		const PS::S64 N_kind[3] = {50, 100, 200}; // ume, shiso, salt
		const PS::S64 N = N_kind[0] + N_kind[1] + N_kind[2];
		ptcl.setNumberOfParticleLocal(N);
		sysinfo.end_time = 1.0;
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
		dinfo.setPosRootDomain(
			PS::F64vec(-kBoxHalfX - 0.002, -kBoxHalfY - 0.002, -kBoxHalfZ - 0.002),
			PS::F64vec( kBoxHalfX + 0.002,  kBoxHalfY + 0.002,  kBoxHalfZ + 0.002)
		);
		if(PS::Comm::getRank() != 0) return;
		const PS::F64 rad_kind[3] = {0.001, 0.0005, 0.00025}; // 1 mm, 0.5 mm, 0.25 mm
		const PS::F64 x_min = -kBoxHalfX;
		const PS::F64 x_max =  kBoxHalfX;
		const PS::F64 y_min = -kBoxHalfY;
		const PS::F64 y_max =  kBoxHalfY;
		const PS::F64 z_min = -kBoxHalfZ;
		const PS::F64 z_max =  kBoxHalfZ;
        const PS::F64 max_occupy_volume = 8*rad_kind[0]*rad_kind[0]*rad_kind[0]*N;
        const PS::F64 ratio = max_occupy_volume / kBoxVolume;
        std::cout<<"max occupy volume ratio: "<<ratio<<std::endl;
        if (ratio*2.0 > 0.99) {
            std::cerr<<"Error: particles are too many to fit in the box."<<std::endl;
            exit(1);
        }
        PS::F64 z_max_init = z_min + (z_max - z_min) * ratio*2.0;
        //if(z_max_init < 0.0) z_max_init = 0.0;
        std::cout<<"z_max= "<<z_max<<" z_max_init: "<<z_max_init<<std::endl;
		int cnt = 0;
		std::vector<PS::F64vec> placed_pos;
		std::vector<PS::F64> placed_rad;
		placed_pos.reserve(N);
		placed_rad.reserve(N);
		for(int k = 0 ; k < 3 ; ++ k){
			const PS::F64 rad = rad_kind[k];
			for(PS::S64 i = 0 ; i < N_kind[k] ; ++ i){
				ThisPtcl ith;
				bool placed = false;
				const int max_trials = 2000;
				for(int trial = 0 ; trial < max_trials && !placed ; ++ trial){
					ith.pos.x = x_min + (x_max - x_min) * rand() / (double)RAND_MAX;
					ith.pos.y = y_min + (y_max - y_min) * rand() / (double)RAND_MAX;
                    ith.pos.z = z_min + (z_max_init - z_min) * rand() / (double)RAND_MAX;
					//PS::F64 z = z_min + (z_max_init - z_min) * rand() / (double)RAND_MAX;
					//ith.pos.z = -std::abs(z); // bias to lower z side
					bool overlap = false;
					for(std::size_t j = 0 ; j < placed_pos.size() ; ++ j){
						const PS::F64vec dp = ith.pos - placed_pos[j];
						const PS::F64 min_dist = rad + placed_rad[j];
						if(dp * dp < min_dist * min_dist){
							overlap = true;
							break;
						}
					}
					if(!overlap){
                        std::cout<<"fail i:"<<i<<std::endl;
						placed = true;
					}
				}
				if(!placed){
					ith.pos.x = 0.0;
					ith.pos.y = 0.0;
					ith.pos.z = z_min + rad;
				}
                //if(ith.pos.z < -0.0045)
                //    std::cout<<"ith.pos= "<<ith.pos<<std::endl;
				ith.id = cnt;
				ith.kind = k;
				ith.rad = rad;
				ith.vel = 0.0;
				ith.avel = 0.0;
				ith.mat  = Test;
				ith.mass = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner = 0.4 * ith.mass * ith.rad * ith.rad;
				ptcl[cnt++] = ith;
				placed_pos.push_back(ith.pos);
				placed_rad.push_back(rad);
			}
		}
	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		const PS::F64vec offset = kVibrationAmp * sin(kVibrationOmega * sysinfo.time);
		const PS::F64vec O = offset;
		InfinitePlane<ThisPtcl> wall_xp(Test, PS::F64vec( 1.0, 0.0, 0.0), O + PS::F64vec( kBoxHalfX, 0.0, 0.0));
		InfinitePlane<ThisPtcl> wall_xm(Test, PS::F64vec(-1.0, 0.0, 0.0), O + PS::F64vec(-kBoxHalfX, 0.0, 0.0));
		InfinitePlane<ThisPtcl> wall_yp(Test, PS::F64vec( 0.0, 1.0, 0.0), O + PS::F64vec(0.0,  kBoxHalfY, 0.0));
		InfinitePlane<ThisPtcl> wall_ym(Test, PS::F64vec( 0.0,-1.0, 0.0), O + PS::F64vec(0.0, -kBoxHalfY, 0.0));
		InfinitePlane<ThisPtcl> wall_zp(Test, PS::F64vec( 0.0, 0.0, 1.0), O + PS::F64vec(0.0, 0.0,  kBoxHalfZ));
		InfinitePlane<ThisPtcl> wall_zm(Test, PS::F64vec( 0.0, 0.0,-1.0), O + PS::F64vec(0.0, 0.0, -kBoxHalfZ));
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			ptcl[i].acc.z -= 9.8;
			wall_xp.setForce(ptcl[i]);
			wall_xm.setForce(ptcl[i]);
			wall_yp.setForce(ptcl[i]);
			wall_ym.setForce(ptcl[i]);
			wall_zp.setForce(ptcl[i]);
			wall_zm.setForce(ptcl[i]);
		}
	}
	static void postTimestep(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
	}
};
