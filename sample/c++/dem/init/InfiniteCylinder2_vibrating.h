#pragma once

#include <vector>

namespace {
	constexpr PS::F64 kCylRadius = 0.005;  // 0.5 cm 円柱の半径
	constexpr PS::F64 kCylLength = 0.05;   // 5 cm 円柱の長さ
	constexpr PS::F64 kYMin = -kCylLength * 0.5;
	constexpr PS::F64 kYMax =  kCylLength * 0.5;
	//constexpr PS::F64vec kVibrationAmp(0.001, 0.0005, 0.001);
	constexpr PS::F64vec kVibrationAmp(0.0, 0.0, kCylRadius*0.1);
	constexpr PS::F64 kVibrationOmega = 2.0 * M_PI; // 1 Hz
}

template <class ThisPtcl> class Problem{
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		const PS::S64 N_kind[3] = {50, 100, 300}; // ume, shiso, salt
		const PS::S64 N = N_kind[0] + N_kind[1] + N_kind[2];
		ptcl.setNumberOfParticleLocal(N);
		//sysinfo.end_time = 10.0;
		sysinfo.end_time = 1.0;
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
		dinfo.setPosRootDomain(PS::F64vec(-0.01, -0.03, -0.01), PS::F64vec(0.01, 0.03, 0.01));
		if(PS::Comm::getRank() != 0) return;
		const PS::F64 rad_kind[3] = {0.001, 0.0005, 0.00025};
		const PS::F64 R = kCylRadius;
		const PS::F64 y_min = kYMin;
		const PS::F64 y_max = kYMax;
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
					const PS::F64 rx = (2.0 * rand() / (double)RAND_MAX - 1.0) * (R - rad);
					PS::F64 rz = (2.0 * rand() / (double)RAND_MAX - 1.0) * (R - rad);
					if(rx * rx + rz * rz > (R - rad) * (R - rad)) continue;
					rz = -std::abs(rz); // bias to lower z side
					ith.pos.x = rx;
					ith.pos.z = rz;
					ith.pos.y = y_min + (y_max - y_min) * rand() / (double)RAND_MAX;
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
						placed = true;
					}
				}
				if(!placed){
					// fallback: accept placement even if overlap to avoid infinite loop
					const PS::F64 rx = 0.0;
					const PS::F64 rz = - (R - rad);
					ith.pos.x = rx;
					ith.pos.z = rz;
					ith.pos.y = y_min + (y_max - y_min) * rand() / (double)RAND_MAX;
				}
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
		const PS::F64vec axis(0.0, 1.0, 0.0);
		const PS::F64 R = kCylRadius;
		const PS::F64vec base_center(0.0, 0.0, 0.0);
		const PS::F64vec base_bottom(0.0, kYMin, 0.0);
		const PS::F64vec offset = kVibrationAmp * sin(kVibrationOmega * sysinfo.time);
		InfiniteCylinder<ThisPtcl> cylinder(Test, axis, base_center + offset, R);
		Disc<ThisPtcl> bottom(Test, axis, base_bottom + offset, 0.0, R);
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			//Earth's gravity
			ptcl[i].acc.z -= 9.8;
			cylinder.setForce(ptcl[i]);
			bottom.setForce(ptcl[i]);
		}
	}
	static void postTimestep(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
	}
};
