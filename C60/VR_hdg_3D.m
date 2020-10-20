%function alphaHull_demo 
% Demo for the alpha hull function 
%
clear; clc;
atom_coor=pdb2mat('c60_Yoshida.pdb'); 
data(:,1)=atom_coor.X;
data(:,2)=atom_coor.Y;
data(:,3)=atom_coor.Z;

[N, ntt]=size(data);
idim=3;

fileID = fopen('./vct_c60rips_eigv.txt','w');



for iloop=1:60

VRdiameter=0.1*iloop;
mat=zeros(N,N);
for i=1:N
    for j=i+1:N
        dist=round(sqrt(sum((data(i,:)-data(j,:)).^2)),2);       
        if(dist<=VRdiameter && i~=j)

            mat(j,i)=1;  
        end
    end
end

[kSkeletonOfVRComplex,simplexDimension]=computeVRComplex(mat,idim);


%%  computer boundaries

bdyMatricies{1}=zeros(N,1);

for it=2:idim+1
    
simplexDimension=it;

KSimplicies=kSkeletonOfVRComplex{simplexDimension};

% if (simplexDimension==1)
%     boundaryMatrix=zeros(size(KSimplicies,1),1)';
%   %  fprintf('simplexdimension: %i %i \n', simplexDimension,it); 
% end

KMinusOneSimplicies=kSkeletonOfVRComplex{simplexDimension-1};
nKSimplicies=size(KSimplicies,1);
nKMinusOneSimplicies=size(KMinusOneSimplicies,1);
boundaryMatrix=zeros(nKMinusOneSimplicies,nKSimplicies);

% fprintf('simplexdimension: %i %i %i\n', simplexDimension,size(boundaryMatrix));
% fprintf('bdmatrix dimension: %i %i %i\n', nKSimplicies,nKMinusOneSimplicies);

for i=1:nKMinusOneSimplicies
   for j=1:nKSimplicies
       
      if KSimplicies(j,1)>KMinusOneSimplicies(i,1)
          break
      end
      
      if KMinusOneSimplicies(i,simplexDimension-1)>KSimplicies(j,simplexDimension)
        continue 
      end
     
       if ismembc(KMinusOneSimplicies(i,:),KSimplicies(j,:))
          [a, indiciesOfMissingElements] = find(ismembc(KSimplicies(j,:), KMinusOneSimplicies(i,:))==0);
          boundaryMatrix(i,j)=(-1)^mod(indiciesOfMissingElements+1,2);  %%%%%%
          %  boundaryMatrix(i,j)=1;  %%%% remove sign
       end
       
   end
end
    bdyMatricies{it}=boundaryMatrix(:,:); 
end

%%   hodge matrix

for id=1:idim
    nmt=size(kSkeletonOfVRComplex{id},1);
    hodgematrix{id}=zeros(nmt,nmt);
    
    matI=bdyMatricies{id};
    matIadd1=bdyMatricies{id+1};
   
    if(id==1)
       hodgematrix{id}=matIadd1*matIadd1';  
    else
       hodgematrix{id}=matI'*matI+matIadd1*matIadd1';
       %newhodgematrix{id}=matI'*matI-matIadd1*matIadd1';  %%%%%%%% L^{low}-L^{upper}
       %newhodgematrix{id}=matIadd1*matIadd1';
       %newhodgematrix{id}=matI'*matI;
    end
end

vct=cell(3,1);
    for i = 1:3
        mathdg=hodgematrix{i};
        %     if numel(mathdg)>0
        [Vec,Dmat] = eig(mathdg);
        ndg=size(Dmat,1);
        ivct=[];
        for indg=1:ndg
            ivct(indg)=Dmat(indg,indg);
            fprintf(fileID,'%3d %4.2f %12.5f\n',i,VRdiameter,ivct(indg));
        end
      end
end
    
fclose(fileID);
