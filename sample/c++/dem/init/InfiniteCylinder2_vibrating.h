#pragma once

template <class ThisPtcl> class Problem{
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		const PS::S64 N = 1000;
		ptcl.setNumberOfParticleLocal(N);
		sysinfo.end_time = 10.0;
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
		dinfo.setPosRootDomain(PS::F64vec(-1.2, -1.2, -1.2), PS::F64vec(1.2, 1.2, 1.2));
		if(PS::Comm::getRank() != 0) return;
		const PS::F64 rad = 0.02;
		const PS::F64 R = 1.0;
		const PS::F64 y_min = -0.8;
		const PS::F64 y_max = 0.8;
		for(PS::S64 i = 0 ; i < N ; ++ i){
			ThisPtcl ith;
			bool placed = false;
			while(!placed){
				const PS::F64 rx = (2.0 * rand() / (double)RAND_MAX - 1.0) * (R - rad);
				const PS::F64 rz = (2.0 * rand() / (double)RAND_MAX - 1.0) * (R - rad);
				if(rx * rx + rz * rz > (R - rad) * (R - rad)) continue;
				ith.pos.x = rx;
				ith.pos.z = rz;
				ith.pos.y = y_min + (y_max - y_min) * rand() / (double)RAND_MAX;
				placed = true;
			}
			ith.id = i;
			ith.rad = rad;
			ith.vel = 0.0;
			ith.avel = 0.0;
			ith.mat  = Test;
			ith.mass = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
			ith.iner = 0.4 * ith.mass * ith.rad * ith.rad;
			ptcl[i] = ith;
		}
	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		const PS::F64vec axis(0.0, 1.0, 0.0);
		const PS::F64 R = 1.0;
		const PS::F64vec base_center(0.0, 0.0, 0.0);
		const PS::F64vec base_bottom(0.0, -1.0, 0.0);
		const PS::F64vec amp(0.1, 0.05, 0.1);
		const PS::F64 omega = 2.0 * M_PI;
		const PS::F64vec offset = amp * sin(omega * sysinfo.time);
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
