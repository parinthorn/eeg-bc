function A=mkfilt_eloreta(L);
% Guido Nolte, 2012-2015
% g.nolte@uke.de

% If you use this code for a publication, please ask Guido Nolte for the
% correct reference to cite.

% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.


[nchan ng ndum]=size(L);
LL=zeros(nchan,ndum,ng);
for i=1:ndum;
    LL(:,i,:)=L(:,:,i);
end
LL=reshape(LL,nchan,ndum*ng);

u0=eye(nchan);
W=reshape(repmat(eye(ndum),1,ng),ndum,ndum,ng);
Winv=zeros(ndum,ndum,ng);
winvkt=zeros(ng*ndum,nchan);
kont=0;
kk=0;
while kont==0;
    kk=kk+1;
    for i=1:ng;
        Winv(:,:,i)=(inv(W(:,:,i)));
        %if i==ng;disp(W(:,:,i));end
    end
    for i=1:ng;
        %winvkt(i,:,:)=Winv(:,:,i)*(squeeze(LL(:,:,i)))';
        %winvkt(i,:,:)=(squeeze(LL(:,:,i)))';
        winvkt(3*(i-1)+1:3*i,:)=Winv(:,:,i)*LL(:,3*(i-1)+1:3*i)';
    end
    kwinvkt=LL*winvkt;
    %kwinvkt(1:4,1:4)
        alpha=.001*trace(kwinvkt)/nchan;
%         alpha=.05*trace(kwinvkt)/nchan;
        M=inv(kwinvkt+alpha*u0);
        
        for i=1:ng;
        Lloc=squeeze(L(:,i,:));
        Wold=W;
        W(:,:,i)=sqrtm(Lloc'*M*Lloc);
        end
    reldef=(norm(reshape(W,[],1)-reshape(Wold,[],1))/norm(reshape(Wold,[],1)));
%     disp(reldef)
    if kk>20 | reldef< .000001 ; kont=1;end;
end
%disp(kk)

ktm=LL'*M;
%ktm=reshape(ktm,ng,ndum,nchan);
 A=zeros(nchan,ng,ndum);

 for i=1:ng;
     A(:,i,:)=(Winv(:,:,i)*ktm(3*(i-1)+1:3*i,:))';
 end
return