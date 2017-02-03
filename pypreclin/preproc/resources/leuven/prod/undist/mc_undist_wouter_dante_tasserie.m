function [] = mc_undist_wouter_dante3(edir,fileName)

dd=dir([edir '*_msk*.nii']);

path_mask=[edir dd(1).name];

V_mask=spm_vol(path_mask);

mask=spm_read_vols(V_mask);

mask=smooth3(mask,'gauss',7);
mask=(mask>0.1*max(mask(:)));

% add black lign with 1 to fit the correct matrix size
if size(mask,1)>size(mask,2)
    mask(:,end+1:end+(size(mask,1)-size(mask,2)),:)= mask(1,1,1);
    V_mask.dim(2)=V_mask.dim(1);
end

% % pos=strfind(path_mask,'_');

dd2=dir([edir fileName]);

path_epi=[edir dd2(1).name];

%path_ave=[path_epi(1:end-4) '_ave.nii'];

V_epi=spm_vol(path_epi);


T=spm_read_vols(V_epi);

if size(T,1)>size(T,2)
    T(:,end+1:end+(size(T,1)-size(T,2)),:,:)= T(1,1,1,1);
    for i=1:size(V_epi)
        V_epi(i).dim(2)=V_epi(i).dim(1);
    end
end

[xdim, ydim zdim tdim]=size(T);


Tn=T;
for k=1:tdim
    vol=squeeze(Tn(:,:,:,k));
    Tn(:,:,:,k)=vol.*mask;
end
Tn(isnan(Tn))=0;  

%T=Tn;

% data rescaling

scaling=1000/max(Tn(:));
T=scaling*T;
clear Tm;


Tave=mean(T,4);

cc=zeros(1,tdim);
for i=1:tdim
    Tc=Tn(:,:,:,i);
    cc(i)=corr(Tc(mask==1),Tave(mask==1));
end
    

[~, nt] = max(cc);


%target_image=T(:,:,:,nt);


seq=100101;
% seq=110111;

options=0;


PE=0;

nx=V_epi(1).dim(1);
ny=V_epi(1).dim(2);
nsl=V_epi(1).dim(3);
ni=length(V_epi);

% % tslc=fix(nsl/2);


for i=1:size(mask,1)
image=mask(i,:,:);
xpixel(i)=sum(image(:));
end

x_list=find(xpixel>max(xpixel)/10);

for i=1:size(mask,2)
image=mask(:,i,:);
ypixel(i)=sum(image(:));
end

y_list=find(ypixel>max(ypixel)/10);

for i=1:size(mask,3)
image=mask(:,:,i);
zpixel(i)=sum(image(:));
end

z_list=find(zpixel>max(zpixel)/20);


marl=x_list(1)-1;

marr=nx-x_list(end);

marb=y_list(1)-1;

mart=ny-y_list(end);

marct=0;

marcb=0;


msg = sprintf('------ Undistortion Parameters ------\n'); disp(msg); 

%msg = sprintf(' Target Dir          = %d\n',tdir); disp(msg); 
msg = sprintf(' UD sequence         = %d\n',seq); disp(msg); 
msg = sprintf(' UD options          = %d\n',options); disp(msg); 
%msg = sprintf(' FFT Scale Factor    = %g\n',fftscale); disp(msg); 

% select slice orientation
%PE=-1 % this is R>>L orientation  ! check !
%PE=0 % this is F>>H orientation  ! check !
%PE=1 % this is L>>R orientation  ! check !
%PE=2 % this is H>>F orientation  ! check !

%
if options==0
% %     doplot1=0 ; % 0-no 1-limited 2-detailed output for individual time points
% %     doplot2=0 ; % summary plot
    dowrite=1 ;
    dooptimize=0;
    doquad=1;
elseif options==1
% %     doplot1=1 ; % 0-no 1-limited 2-detailed output for individual time points
% %     doplot2=1 ; % summary plot
    dowrite=1 ;
    dooptimize=1;
    doquad=1;
elseif options==2 % same as options=0 but 4 quadrants
% %     doplot1=0 ; % 0-no 1-limited 2-detailed output for individual time points
% %     doplot2=0 ; % summary plot
    dowrite=1 ;
    dooptimize=0;
    doquad=4;
elseif options==3 % same as options=1 but 4 quadrants
% %     doplot1=1 ; % 0-no 1-limited 2-detailed output for individual time points
% %     doplot2=1 ; % summary plot
    dowrite=1 ;
    dooptimize=1;
    doquad=4;
end
%
% if dooptimize==1
%     slin=tslc;
%     slst=1;
%     slfi=slin;
% else
%     slin=1;
%     slst=1;
%     slfi=nsl;
% end


%%%

%msg = sprintf('\n### Starting Undistortion Correction ###'); disp(msg);    


%%% INFO FOR HEADER FILE
% Nx=nx;
% Ny=ny;
% ni=Rep;  
% rnx=nx;
% rny=ny;
% 
% ResX = FovX/Nx;
% ResY = FovY/Ny;
% ResZ = dSlice;

%%%%%%%%%%
% Read target volume



%%%%%%%%%%


    %Ttemp=R(31:51,31:51,25)./T(31:51,31:51,25,1);
    %T=T.*mean(Ttemp(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over nr images in run
    nr=ni;
    %
% %     diff1=zeros(nr,1);
% %     diff2=zeros(nr,1);
% %     diff3=zeros(nr,1);
    %
    %M=zeros(nx,ny,nsl,ni);
    M=T;
    
    if PE==2||PE==0
        % temp output array
        rnx=nx;
        rny=ny;
    elseif PE==1||PE==-1
        rnx=ny;
        rny=nx;
    end

    out=zeros(nx,ny,nr);

    
    % smoothing

    psf=fspecial('gaussian',5,1);

    % init statistics vars

    %%%%%%%%%%
    % slice loop start

    for w=1:length(z_list) %  1:nsl
        
        slc=z_list(w);
        
        msg = sprintf('undist. slice : %d',slc); disp(msg);    

        
        %%%%%%%%%%
        % time loop start

        for nc=1:ni
        %for nc=25:34
            % align data images 
            
           
        ttmp=squeeze(T(:,:,slc,nt));
        
        lim=max(ttmp(:));
        clim=[0 lim]*3/2;
% %         dlim=[-lim lim];

        if PE==2
            im1=ttmp;
        elseif PE==-1
            im1=rot90( ttmp ,1 );
        elseif PE==0
            im1=rot90( ttmp ,2 );
        elseif PE==1
            im1=rot90( ttmp ,3 );
        end

               
% %         % smoothed target image
% %         im=imfilter(im1,psf,'same','conv');
% %         tim1=deconvlucy(edgetaper(im,psf),psf,50);    
        tim1=imfilter(im1,psf,'same','conv', 'replicate');

        % determine avg. image noise level and threshold
        taxi=rnx-10;
        taxf=rnx;
        tayi=rny-10;
        tayf=rny;
        tltmp=abs(im1(taxi:taxf,tayi:tayf));
        tlev=sum(tltmp(:))/10/10*1.25;
        
        % create mask1
        mask1=ones(size(im1));
        mask1(im1(:)<tlev)=0;
%        mask1(:,1:20)=0;
        
        tim1=tim1.*mask1;

        %
        qnn=doquad;
        stb00=zeros(ni,qnn);
        stb11=zeros(ni,qnn);
% %         stb21=zeros(ni,qnn);
% %         stb02=zeros(ni,qnn);
% %         stbx2=zeros(ni,qnn);
% %         stby2=zeros(ni,qnn);
        stb02=zeros(ni,qnn);
% %         stco=zeros(ni,qnn);
% %         stn=[1:ni];

                       
            
            
            tmp=T(:,:,slc,nc);
% %             tmp=Tx(:,:,slc,nc);
            if PE==2
                im2=tmp;
            elseif PE==-1rnx
                im2=rot90( tmp ,1 );
            elseif PE==0
                im2=rot90( tmp ,2 );
            elseif PE==1
                im2=rot90( tmp ,3 );
            end

            % smoothed distorted image
% %             im=imfilter(im2,psf,'same','conv');
% %             tim2=deconvlucy(edgetaper(im,psf),psf,20);
            tim2=imfilter(im2,psf,'same','conv', 'replicate');
% %             tim2=im2;

            % create mask2
            mask2=ones(size(im2));
            mask2(tim2(:)<tlev)=0;
%            mask2(:,1:20)=0;

            % distorted image
            tim2=tim2.*mask2;
                
            % init final image
            im3=zeros(size(im2));

            %
            tmask = mask1;

            % init structure for distortion correction
            V.tmask=tmask;
            V.nx=rnx;
            V.ny=rny;
            W.bim=zeros(size(im1));

            qim1=zeros(size(tim1));
            qim1(:,:)=tim1(:,:);
            qim2=zeros(size(tim2));
            qim2(:,:)=tim2(:,:);

            Wqim=zeros(size(im1));

            if doquad==4
                qim1=zeros(size(tim1,1),size(tim1,2),4);
                qim1(1:end/2,1:end/2,1)=tim1(1:end/2,1:end/2);
                qim1(1:end/2,end/2+1:end,2)=tim1(1:end/2,end/2+1:end);
                qim1(end/2+1:end,end/2+1:end,3)=tim1(end/2+1:end,end/2+1:end);
                qim1(end/2+1:end,1:end/2,4)=tim1(end/2+1:end,1:end/2);
                qim2=zeros(size(tim2,1),size(tim2,2),4);
                qim2(1:end/2,1:end/2,1)=tim2(1:end/2,1:end/2);
                qim2(1:end/2,end/2+1:end,2)=tim2(1:end/2,end/2+1:end);
                qim2(end/2+1:end,end/2+1:end,3)=tim2(end/2+1:end,end/2+1:end);
                qim2(end/2+1:end,1:end/2,4)=tim2(end/2+1:end,1:end/2);
                qim3=zeros(size(tim1));
            end

            for qq=1:qnn
            
                if (length(find(tmask))>100);


                    tqim1=qim1(:,:,qq);
                    tqim2=qim2(:,:,qq);
                    
                    % target image
                    V.im1=zeros(size(tqim1));
    %                V.im1(marl+1:end-marr,marb+1:end-mart)=tim1(marl+1:end-marr,marb+1:end-mart).*tim1(marl+1:end-marr,marb+1:end-mart);
                    V.im1(marl+1:end-marr,mart+1:end/2-marct)=tqim1(marl+1:end-marr,mart+1:end/2-marct).*tim1(marl+1:end-marr,mart+1:end/2-marct);
                    V.im1(marl+1:end-marr,end/2+marcb+1:end-marb)=tqim1(marl+1:end-marr,end/2+marcb+1:end-marb).*tim1(marl+1:end-marr,end/2+marcb+1:end-marb);


                    % distorted image
                    V.im2=zeros(size(tqim2));
    %                V.im2(marl+1:end-marr,marb+1:end-mart)=tim2(marl+1:end-marr,marb+1:end-mart).*tim2(marl+1:end-marr,marb+1:end-mart);
                    V.im2(marl+1:end-marr,mart+1:end/2-marct)=tqim2(marl+1:end-marr,mart+1:end/2-marct).*tim2(marl+1:end-marr,mart+1:end/2-marct);
                    V.im2(marl+1:end-marr,end/2+marcb+1:end-marb)=tqim2(marl+1:end-marr,end/2+marcb+1:end-marb).*tim2(marl+1:end-marr,end/2+marcb+1:end-marb);


    %                subplot(3,2,1);
    %                imagesc(abs(tim1'));
    %                subplot(3,2,2);
    %                imagesc(abs(V.im1'));
    %                drawnow;

                    %%% vars for optimization loop

                    V.scale_list = [0 64 32 16 8 4 2 1];

                    % order in lists: [b00 b01 b02 b10 b11 b20]
                    % b00: shift in PE
                    % b01: 1st order modulation in RO of shift in PE 
                    % b02: 2nd order modulation in RO of shift in PE
                    % b10: 1st order distortion in PE
                    % b11: 1st order modulation in RO of distortion in PE
                    % b20: 2nd order distortion in PE

                    V.db.step = 0.010*[1 1 1 1 1 1];
                    %  V.db.step = [0.025 0.025 0.025 0.025 0.025 0.025];

                    % define type of correction here 
                    seq_list=[0 0 0 0 0 0];
                    tst=seq;
                    for ii=0:5
                        if (tst-10^(5-ii))>=0
                            seq_list(ii+1)=1;
                            tst=tst-10^(5-ii);
                        end
                    end

                    V.db.dolist = seq_list;

    %                if (PE==0) || (PE==2) % F-H or H-F
                        % focus on distortion
    %                    V.db.dolist = [1 0 0 1 0 1];
    %                    if (slc>30)
    %                        V.db.dolist = [1 0 0 1 0 0];
    %                    end
                        %
    %                elseif (PE==-1) || (PE==1) % R-L or L-R
                        % focus on shifts
    %                    V.db.dolist = [1 1 1 1 0 0];
                        %
    %                end

                    V.db.list1 = 0; if V.db.dolist(1); V.db.list1=[-1 0 1]; end
                    V.db.list2 = 0; if V.db.dolist(2); V.db.list2=[-1 0 1]; end
                    V.db.list3 = 0; if V.db.dolist(3); V.db.list3=[-1 0 1]; end
                    V.db.list4 = 0; if V.db.dolist(4); V.db.list4=[-1 0 1]; end
                    V.db.list5 = 0; if V.db.dolist(5); V.db.list5=[-1 0 1]; end
                    V.db.list6 = 0; if V.db.dolist(6); V.db.list6=[-1 0 1]; end


    %%                  % subroutine

    %                    tic;
% %                         W=mc_undist_loop_wouter(V) ;
                        W=mc_undist_loop_wouter_fast_tasserie(V) ; % updated (08/06/2011)
    %                    toc;



                    % ouput
                    % im3=W.im3;
                    stb00(nc,qq)=W.b00;
                    stb01(nc,qq)=W.b01;
                    stb02(nc,qq)=W.b02;
                    stb10(nc,qq)=W.b10;
                    stb11(nc,qq)=W.b11;
                    stb20(nc,qq)=W.b20;
% %                     stcst(nc,qq)=W.cst;
                    scale=W.scale;
                    iter=W.iter;
                    
                    if options==1;
                    fprintf(1,'rep = %d: scale = %d (%d): b00 = %1.3f, b01 = %1.3f, b02 = %1.3f, b10 = %1.3f, b11 = %1.3f, b20 = %1.3f, cost = %1.5f\n',nc,scale,iter,W.b00,W.b01,W.b02,W.b10,W.b11,W.b20,W.cst*1e16);
                    end

                    [X,Y] = meshgrid(1:V.ny,1:V.nx);
                    qim3 = interp2(im2,X+W.bim,Y,'cubic');

                else % tmask

% %                     fprintf(1,'rep = %d \n',nc);
                    qim3=im1;
                
                end % tmask

                if doquad==1
                    im3=qim3;
                elseif doquad==4
% %                 pn=0;
                px=2;
                py=3;
                pn=1;
                subplot(py,px,pn);
                if qq==1
                        im3(1:end/2,1:end/2)=qim3(1:end/2,1:end/2);
                        Wqim(1:end/2,1:end/2)=W.bim(1:end/2,1:end/2);
                    elseif qq==2
                        im3(1:end/2,end/2+1:end)=qim3(1:end/2,end/2+1:end);
                        Wqim(1:end/2,end/2+1:end)=W.bim(1:end/2,end/2+1:end);
                    elseif qq==3
                        im3(end/2+1:end,end/2+1:end)=qim3(end/2+1:end,end/2+1:end);
                        Wqim(end/2+1:end,end/2+1:end)=W.bim(end/2+1:end,end/2+1:end);
                    elseif qq==4
                        im3(end/2+1:end,1:end/2)=qim3(end/2+1:end,1:end/2);
                        Wqim(end/2+1:end,1:end/2)=W.bim(end/2+1:end,1:end/2);
                    end
                end            
            end % qq
            
            im3(isnan(im3))=0;
            
            %%%  rotate i=mage back to original direction
            if PE==2
            out(:,:,nc)=im3;
            elseif PE==-1
            out(:,:,nc)=rot90( im3 ,3);
            elseif PE==0
            out(:,:,nc)=rot90( im3 ,2);
            elseif PE==1
            out(:,:,nc)=rot90( im3 ,1 );
            end

            M(:,:,slc,nc)=out(:,:,nc);

            if dooptimize==1
                colormap(jet);
% %                 pn=0;
                px=2;
                py=3;
                %
%                pn=pn+1;
%                subplot(py,px,pn);
%                imagesc(abs(im1'));

                %
                pn=1;
                subplot(py,px,pn);
                imagesc(Wqim',[-3 3]);
                %
                pn=2;
                subplot(py,px,pn);
                llm=2; ulm=llm*3;
                if qnn==1
                    xl=1:nc;
                    plot(xl,stb00(xl),'.',xl,stb01(xl),'.',xl,stb02(xl),'.',xl,stb10(xl),'.',xl,stb11(xl),'.',xl,stb20(xl),'.')
                    set(gca,'xLim',[0 ni+1]);
                    set(gca,'yLim',[-llm llm]);
                elseif qnn==4
                    xl=1:nc;
                    plot( ...
                        xl,stb00(xl,1)+ulm,'.',xl,stb01(xl,1)+ulm,'.',xl,stb02(xl,1)+ulm,'.',xl,stb10(xl,1)+ulm,'.',xl,stb11(xl,1)+ulm,'.',xl,stb20(xl,1)+ulm,'.', ...
                        xl,stb00(xl,2)+llm,'.',xl,stb01(xl,2)+llm,'.',xl,stb02(xl,2)+llm,'.',xl,stb10(xl,2)+llm,'.',xl,stb11(xl,2)+llm,'.',xl,stb20(xl,2)+llm,'.', ...
                        xl,stb00(xl,3)-llm,'.',xl,stb01(xl,3)-llm,'.',xl,stb02(xl,3)-llm,'.',xl,stb10(xl,3)-llm,'.',xl,stb11(xl,3)-llm,'.',xl,stb20(xl,3)-llm,'.', ...
                        xl,stb00(xl,4)-ulm,'.',xl,stb01(xl,4)-ulm,'.',xl,stb02(xl,4)-ulm,'.',xl,stb10(xl,4)-ulm,'.',xl,stb11(xl,4)-ulm,'.',xl,stb20(xl,4)-ulm,'.' ...
                        );
                    set(gca,'xLim',[0 ni+1]);
                    set(gca,'yLim',[-(llm+ulm) llm+ulm]);
                end
                %
            	xl=nx+1;
                yl=ny+1;
                %
                pn=pn+1;
                subplot(py,px,pn);
                imagesc(abs(tim2'),clim*2/3);
                line([marl marl],[0 yl],'Color',[0 0 0],'Linewidth',2)
                line([rnx-marr rnx-marr],[0 yl],'Color',[0 0 0],'Linewidth',2)
                line([0 xl],[mart mart],'Color',[0 0 0],'Linewidth',2)
                line([0 xl],[rny-marb rny-marb ],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2+marcb rny/2+marcb],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2-marct rny/2-marct ],'Color',[0 0 0],'Linewidth',2)
                %
                pn=pn+1;
                subplot(py,px,pn);
                imagesc(abs(im3'),clim*2/3);
                line([marl marl],[0 yl],'Color',[0 0 0],'Linewidth',2)
                line([rnx-marr rnx-marr],[0 yl],'Color',[0 0 0],'Linewidth',2)
                line([0 xl],[mart mart],'Color',[0 0 0],'Linewidth',2)
                line([0 xl],[rny-marb rny-marb ],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2+marcb rny/2+marcb],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2-marct rny/2-marct ],'Color',[0 0 0],'Linewidth',2)
                %
                pn=pn+1;
                subplot(py,px,pn);
                imagesc(abs(tim1-tim2)',clim/5);
                line([marl marl],[0 yl],'Color',[0 0 0],'Linewidth',2);
                line([rnx-marr rnx-marr],[0 yl],'Color',[0 0 0],'Linewidth',2);
                line([0 xl],[mart mart],'Color',[0 0 0],'Linewidth',2)
                line([0 xl],[rny-marb rny-marb ],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2+marcb rny/2+marcb],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2-marct rny/2-marct ],'Color',[0 0 0],'Linewidth',2)
                %
                pn=pn+1;
                subplot(py,px,pn);
                imagesc(abs(im1-im3)',clim/5);
                line([marl marl],[0 yl],'Color',[0 0 0],'Linewidth',2);
                line([rnx-marr rnx-marr],[0 yl],'Color',[0 0 0],'Linewidth',2);
                line([0 xl],[mart mart],'Color',[0 0 0],'Linewidth',2)
                line([0 xl],[rny-marb rny-marb ],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2+marcb rny/2+marcb],'Color',[0 0 0],'Linewidth',2)
%                line([0 xl],[rny/2-marct rny/2-marct ],'Color',[0 0 0],'Linewidth',2)
                %
%            pause
            drawnow;
            end

        end % nc

    end % slc
    


    if dowrite==1

        
        path_und=[edir 'u_' dd2(1).name];
        
        for i=1:length(V_epi)
        V_epi(i).fname=path_und;
        spm_write_vol(V_epi(i),squeeze(M(:,:,:,i)));
        end

    end % dowrite
    
    
 

% statistics on noise reduction
target_image_pre=squeeze(T(:,:,:,1));
target_image_post=squeeze(M(:,:,:,1));

corr_pre=zeros(1,tdim);
corr_post=zeros(1,tdim);
for k=1:tdim
    im_pre=squeeze(T(:,:,:,k));
    im_post=squeeze(M(:,:,:,k));
    corr_pre(k)=corr(target_image_pre(mask==1),im_pre(mask==1));
    corr_post(k)=corr(target_image_post(mask==1),im_post(mask==1));    
end

figure('visible','off'); plot(corr_pre); hold on; plot(corr_post,'r'); legend('before correction','after correction');
print([path_und(1:end-4) '_stats.tif'],'-dtiff','-r150');

save([path_und(1:end-4) '_stats.mat'],'corr_pre','corr_post');
