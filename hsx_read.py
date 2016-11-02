# Read HSX content from Siesta *.HSX file
# Written by David Visontai, 2016

import HSX_Mod

RytoeV = = 1./13.60580	# Rydberg to electron Volt 

def read_hsx(hsx_filename):
	hsx = HSX_Mod.read_hsx_file(hsx_filename)	
	return hsx


systemlabel="emol"
## READ FILE "system_label.HSX"
hsx_filename =  systemlabel + ".HSX"
hsx = read_hsx(hsx_filename)

#Set units to eV
hsx.hamilt = hsx.hamilt/RytoeV;

gspin = hsx.nspin
if hsx.nspin > 2:
	gspin = 1

# CREATE ORBITAL INFO
iorb = []
#for iuo=1:hsx.no_u
#	iat    = hsx.iaorb(iuo);
#  	iphorb = hsx.iphorb(iuo);
#  	spec   = hsx.isa(iat);
#  	elmti  = elmtxv(iat);
#  	dum    = [ iat, hsx.nquant(spec,iphorb), hsx.lquant(spec,iphorb),0,hsx.zeta(spec,iphorb), coords(iat,1:3), elmti, atmass(elmti)];
#  	inputP.iorb = [inputP.iorb; dum];

#if hsx.nspin > 2:
#	hsx.iorb = [inputP.iorb; inputP.iorb]

tt=[]; tk=[];
for iuo=1:inputP.no_u
	tt = [ tt 1:hsx.numh(iuo)];
	tk= [tk hsx.listhptr(iuo).*ones(hsx.numh(iuo),1)'];%'
end
ind=(tt+tk)';%'

iuo=[];
for ii=1:inputP.no_u
	iuo = [iuo ii.*ones(hsx.numh(ii),1)'];%'
end
juo = (hsx.indxuo(hsx.listh(ind')))';

%% COMPUTE H and S FOR EXTENDED MOLECULE
Rbloch = hsx.xij(ind,1:3) + coords(inputP.iorb(iuo',1),1:3) - coords(inputP.iorb(juo',1),1:3);
HS=[];
for ik=1:size(inputP.EM.kpoints,1)
	eikr = exp(-1.0i* (inputP.EM.kpoints(ik,1).*Rbloch(ind,1)+inputP.EM.kpoints(ik,2).*Rbloch(ind,2)));
  	if inputP.nspin <= 2
  		HS(ik,1).S = sparse(juo,iuo, hsx.Sover.*eikr,inputP.no_u,inputP.no_u);
		for is=1:inputP.nspin
 	    	HS(ik,is).H = sparse(juo,iuo, hsx.hamilt(:,is).*eikr,inputP.no_u,inputP.no_u);
	   	end
    elseif inputP.nspin == 4
    	HS(ik,1).S = kron(eye(2),sparse(juo,iuo, hsx.Sover.*eikr,inputP.no_u,inputP.no_u));
 	    H1=kron([1  0; 0 0], sparse(juo,iuo, hsx.hamilt(:,1).*eikr,inputP.no_u,inputP.no_u));
 	    H2=kron([0  0; 0 1], sparse(juo,iuo, hsx.hamilt(:,2).*eikr,inputP.no_u,inputP.no_u));
 	    H3=kron([0  1; 0 0], sparse(juo,iuo, (hsx.hamilt(:,3)-1i*hsx.hamilt(:,4)).*eikr,inputP.no_u,inputP.no_u));
 	    H4=kron([0  0; 1 0], sparse(juo,iuo, (hsx.hamilt(:,3)+1i*hsx.hamilt(:,4)).*eikr,inputP.no_u,inputP.no_u));
 	    HS(ik,1).H = H1 + H2 + H3 + H4; 
 	elseif inputP.nspin == 8
    	HS(ik,1).S = kron(eye(2),sparse(juo,iuo, hsx.Sover.*eikr,inputP.no_u,inputP.no_u));
 	    H1=kron([1  0; 0 0], sparse(juo,iuo, (hsx.hamilt(:,1)+1i*hsx.hamilt(:,5)).*eikr,inputP.no_u,inputP.no_u));
 	    H2=kron([0  0; 0 1], sparse(juo,iuo, (hsx.hamilt(:,2)+1i*hsx.hamilt(:,6)).*eikr,inputP.no_u,inputP.no_u));
  	    H3=kron([0  1; 0 0], sparse(juo,iuo, (hsx.hamilt(:,3)-1i*hsx.hamilt(:,4)).*eikr,inputP.no_u,inputP.no_u));
 	    H4=kron([0  0; 1 0], sparse(juo,iuo, (hsx.hamilt(:,7)+1i*hsx.hamilt(:,8)).*eikr,inputP.no_u,inputP.no_u));
 	    HS(ik,1).H = H1 + H2 + H3 + H4; 
    end
end

end 
