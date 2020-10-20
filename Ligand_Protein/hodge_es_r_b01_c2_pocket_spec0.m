listfile=fopen('2016_new_2.txt');
indexlist=textscan(listfile,'%s');
for iname =1:100
    indexname=indexlist{1,1}{iname};
    data = pdb2mat1(strcat(indexname,'_pocket.pdb'));
    fileID= fopen(strcat(indexname,'_ligands.mol2'));
    C=textscan(fileID,'%d %s %f %f %f %s %d %s %f','HeaderLines',1,'Delimiter',' ','MultipleDelimsAsOne',1);
    Protein_ele=["C","N","O","S"];
    Ligand_ele=["C","N","P","O","S","F","Cl","Br","I"];
    idim=2;
    nloop=250;
    file0 = fopen(strcat('./',indexname,'_ES_b01_c2_pocket_eigv.txt'),'w');
    file1 = fopen(strcat('./',indexname,'_ES_b01_c2_pocket_simp.txt'),'w');
    for eleIndex_P=1:length(Protein_ele)
        for eleIndex_L=1:length(Ligand_ele)


            b = 1;
            for a=1:length(data.atomName)
                if strcmp(data.atomName{a}(1),Protein_ele(eleIndex_P))
                    atom(b,1)=data.X(a);
                    atom(b,2)=data.Y(a);
                    atom(b,3)=data.Z(a);
                    b=b+1;
                end
            end
            NPro = b-1;
            b=1;
            for a=1:length(C{1,3})
                truncatedLEle=split(C{1,6}(a),'.');
                if strcmp(truncatedLEle{1},Ligand_ele(eleIndex_L))
                    atom(NPro+b,1)=C{1,3}(a);
                    atom(NPro+b,2)=C{1,4}(a);
                    atom(NPro+b,3)=C{1,5}(a);
                    b=b+1;
                end
            end
            NLig=b-1;

            N= NPro + NLig;
            if NPro ==0 || NLig==0 || N<4
                continue
            end
           %% main part, build Rips complex based on modified M matrix
            eleIndex=(eleIndex_P-1)*9+eleIndex_L;
            for iloop=1:nloop
                VRdiameter=0.1*iloop;
                mat=zeros(N,N);
                for i=1:NPro
                    for j=NPro+1:N
                        dist=sqrt(sum((atom(i,:)-atom(j,:)).^2));       
                        if(dist<VRdiameter)

                            mat(j,i)=1;  
                        end
                    end
                end

                [kSkeletonOfVRComplex,simplexDimension]=computeVRComplex(mat,idim);

               %%  computer boundaries

                bdyMatricies{1}=zeros(N,1);
                for it=2:idim     
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
               %% calculate eigenvectors
                for id=1:idim-1
                    nmt=size(kSkeletonOfVRComplex{id},1);
                    fprintf(file1,'%5d %5d %6.2f %7d\n',eleIndex,id,VRdiameter,nmt);
                    hodgematrix{id}=zeros(nmt,nmt);
                    matI=bdyMatricies{id};
                    matIadd1=bdyMatricies{id+1};
                    if(id==1)
                        hodgematrix{id}=matIadd1*matIadd1';
                    else
                        hodgematrix{id}=matI'*matI+matIadd1*matIadd1';
                    end
                end
                for i =1:idim-1
                    mathdg=hodgematrix{i};
                %     if numel(mathdg)>0
                    [Vec,Dmat] = eig(mathdg);
                    ndg=size(Dmat,1);
                    ivct=[];
                    for indg=1:ndg
                        ivct(indg)=Dmat(indg,indg);
                        fprintf(file0,'%5d %5d %6.2f %18.3e\n',eleIndex,i,VRdiameter,ivct(indg));
                    end
                end

            end

 
            clearvars -except indexlist fileID C idim data Protein_ele Ligand_ele eleIndex_P eleIndex_L indexname nloop file0 file1;

        end

    end
    fclose(file0);
    fclose(file1);
end

