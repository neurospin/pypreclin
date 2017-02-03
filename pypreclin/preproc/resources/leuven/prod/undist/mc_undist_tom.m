function [] = mc_undist(ddir)

%% Read parameters from files 

fid = fopen('./prod.par','rb');

% catch problems with old par files
seq=0;
fftscale=1.5e9;
TE = 19;
TR = 2000;
tslc=25;
marl= 0;
marr= 0;
mart= 0;
marb= 0;
marct= 0;
marcb= 0;
%
sLine=fgets(fid);
count=0;
while (sLine>-1)
    sLine=fgets(fid);

    eq=findstr(sLine,'=');
    if(~isempty(eq))
        sParameter=deblank(sLine(1:eq-1));
        psRemainder=sLine(eq+1:end);
        count=count+1;

        %//
        %// Pick the values of important parameters.
        %//

        %data        
        if ( strcmp(sParameter,'Repetitions') )
             Rep = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'Slices') )
             nsl = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'PE direction EPI') )
             PE = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'Matrix X') )
             nx = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'FOV X') )
             FovX = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'Matrix Y') )
             ny = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'FOV Y') )
             FovY = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'Slice Thickness') )
             dSlice = sscanf(psRemainder,'%f');
        elseif ( strcmp(sParameter,'TE') )
             TE = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'TR') )
             TR = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'Target directory') )
             tdir = sscanf(psRemainder,'%s');
        elseif( strcmp(sParameter,'Target image #') )
             nt = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'Target slice #') )
             tslc = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'Margin left') )
             marl = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'Margin right') )
             marr = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'Margin top') )
             mart = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'Margin bottom') )
             marb = sscanf(psRemainder,'%d');
        elseif ( strcmp(sParameter,'Margin center-top') )
             marct = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'Margin center-bot') )
             marcb = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'UD sequence') )
             seq = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'UDoptions') )
             options = sscanf(psRemainder,'%d');
        elseif( strcmp(sParameter,'FFT Scale Factor') )
             fftscale = sscanf(psRemainder,'%f');
        end
    end
end
status=fclose(fid);

msg = sprintf('------ Undistortion Parameters ------\n'); disp(msg); 

msg = sprintf(' Target Dir          = %d\n',tdir); disp(msg); 
msg = sprintf(' UD sequence         = %d\n',seq); disp(msg); 
msg = sprintf(' UD options          = %d\n',options); disp(msg); 
msg = sprintf(' FFT Scale Factor    = %g\n',fftscale); disp(msg); 

% select slice orientation
%PE=-1 % this is R>>L orientation  ! check !
%PE=0 % this is F>>H orientation  ! check !
%PE=1 % this is L>>R orientation  ! check !
%PE=2 % this is H>>F orientation  ! check !

%
if options==0
    doplot1=0 ; % 0-no 1-limited 2-detailed output for individual time points
    doplot2=0 ; % summary plot
    dowrite=1 ;
    dooptimize=0;
elseif options==1
    doplot1=1 ; % 0-no 1-limited 2-detailed output for individual time points
    doplot2=1 ; % summary plot
    dowrite=0 ;
    dooptimize=1;
end
%
if dooptimize==1
    slin=tslc;
    slst=1;
    slfi=slin;
else
    slin=1;
    slst=1;
    slfi=nsl;
end


%%%

%msg = sprintf('\n### Starting Undistortion Correction ###'); disp(msg);    


%%% INFO FOR HEADER FILE
Nx=nx;
Ny=ny;
ni=Rep;  
rnx=nx;
rny=ny;

ResX = FovX/Nx;
ResY = FovY/Ny;
ResZ = dSlice;

%%%%%%%%%%
% Read target volume

T=zeros(nx,ny,nsl,ni);
msg = sprintf('\n- Reading target data'); disp(msg);    
tfile =strcat(tdir,'f.bfloat')
tfid=fopen(tfile,'rb');
tf = fread(tfid,nx*ny*nsl*ni,'float32');
size(tf);
nx*ny*nsl*ni;
status=fclose(tfid);

T = reshape(tf,nx,ny,nsl,ni);
%T=T*(4000/max(max(max(T))));
%%%%%%%%%%
% loop over run directories

    %
    % Read image data
    R=zeros(nx,ny,nsl,ni);

    file =strcat('f.bfloat')
    fid=fopen(file,'rb');
    f = fread(fid,nx*ny*nsl*ni,'float32');
    size(f);
    status=fclose(fid);
    R = reshape(f,nx,ny,nsl,ni); 
    %R=R*(4000/max(max(max(R))));
  
    %Ttemp=R(31:51,31:51,25)./T(31:51,31:51,25,1);
    %T=T.*mean(Ttemp(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over nr images in run
    nr=ni;
    %
    diff1=zeros(nr,1);
    diff2=zeros(nr,1);
    diff3=zeros(nr,1);
    %
    M=zeros(nx,ny,nsl,ni);

    if PE==2|PE==0
        % temp output array
        rnx=nx;
        rny=ny;
    elseif PE==1|PE==-1
        rnx=ny;
        rny=nx
    end

    out=zeros(nx,ny,nr);

    
    % smoothing

    psf=fspecial('gaussian',5,1);

    % init statistics vars

    %%%%%%%%%%
    % slice loop start

    for slc=slin:slst:slfi

        msg = sprintf('undist. slice : %d',slc); disp(msg);    

        stb00=zeros(1,ni);
        stb11=zeros(1,ni);
        stb21=zeros(1,ni);
        stb02=zeros(1,ni);
        stbx2=zeros(1,ni);
        stby2=zeros(1,ni);
        stb02=zeros(1,ni);
        stco=zeros(1,ni);
        stn=[1:ni];

        ttmp=squeeze(T(:,:,slc,nt));

        
        %lim=4000;
        lim=max(ttmp(:));
        clim=[0 lim]*3/2;
        dlim=[-lim lim];

        if PE==2
            im1=ttmp;
        elseif PE==-1
            im1=rot90( ttmp ,1 );
        elseif PE==0
            im1=rot90( ttmp ,2 );
        elseif PE==1
            im1=rot90( ttmp ,3 );
        end

               
        % smoothed target image
        im=imfilter(im1,psf,'same','conv');
        tim1=deconvlucy(edgetaper(im,psf),psf,50);    

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
%       mask1(:,1:20)=0;
        
        tim1=tim1.*mask1;
        
        %%%%%%%%%%
        % time loop start

        for nc=1:ni
        %for nc=25:34
            % align data images 
            tmp=R(:,:,slc,nc);
            if PE==2
                im2=tmp;
            elseif PE==-1
                im2=rot90( tmp ,1 );
            elseif PE==0
                im2=rot90( tmp ,2 );
            elseif PE==1
                im2=rot90( tmp ,3 );
            end

            % smoothed distorted image
            im=imfilter(im2,psf,'same','conv');
            tim2=deconvlucy(edgetaper(im,psf),psf,20);

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
            W.bim=size(nx,ny);

            if (length(find(tmask))>100);

                % target image
                V.im1=zeros(size(tim1));
%                V.im1(marl+1:end-marr,marb+1:end-mart)=tim1(marl+1:end-marr,marb+1:end-mart).*tim1(marl+1:end-marr,marb+1:end-mart);
                V.im1(marl+1:end-marr,mart+1:end/2-marct)=tim1(marl+1:end-marr,mart+1:end/2-marct).*tim1(marl+1:end-marr,mart+1:end/2-marct);
                V.im1(marl+1:end-marr,end/2+marcb+1:end-marb)=tim1(marl+1:end-marr,end/2+marcb+1:end-marb).*tim1(marl+1:end-marr,end/2+marcb+1:end-marb);

                
                % distorted image
                V.im2=zeros(size(tim2));
%                V.im2(marl+1:end-marr,marb+1:end-mart)=tim2(marl+1:end-marr,marb+1:end-mart).*tim2(marl+1:end-marr,marb+1:end-mart);
                V.im2(marl+1:end-marr,mart+1:end/2-marct)=tim2(marl+1:end-marr,mart+1:end/2-marct).*tim2(marl+1:end-marr,mart+1:end/2-marct);
                V.im2(marl+1:end-marr,end/2+marcb+1:end-marb)=tim2(marl+1:end-marr,end/2+marcb+1:end-marb).*tim2(marl+1:end-marr,end/2+marcb+1:end-marb);

                
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


                %% subroutine start here 
                
                % !! WORK IN PROGRESS !! %
%                tic;
                W=mc_undist_loop(V) ;
%                toc;
                % ouput 
                % im3=W.im3;
                stb00(nc)=W.b00;
                stb01(nc)=W.b01;
                stb02(nc)=W.b02;
                stb10(nc)=W.b10;
                stb11(nc)=W.b11;
                stb20(nc)=W.b20;
                stcst(nc)=W.cst;
                scale=W.scale;
                iter=W.iter;

                if options==1;
                fprintf(1,'rep = %d: scale = %d (%d): b00 = %1.3f, b01 = %1.3f, b02 = %1.3f, b10 = %1.3f, b11 = %1.3f, b20 = %1.3f, cost = %1.5f\n',nc,scale,iter,W.b00,W.b01,W.b02,W.b10,W.b11,W.b20,W.cst*1e16);
                end
                
                [X,Y] = meshgrid(1:V.ny,1:V.nx);
                im3 = interp2(im2,X+W.bim,Y,'cubic');

            else % tmask

                fprintf(1,'rep = %d \n',nc);
                im3=im1;

            end % tmask

            im3(find(isnan(im3)))=0;
    
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
                pn=0;
                px=2;
                py=3;
                %
%                pn=pn+1;
%                subplot(py,px,pn);
%                imagesc(abs(im1'));

                %
                pn=pn+1;
                subplot(py,px,pn);
                imagesc(W.bim',[-3 3]);
                %
                pn=pn+1;
                subplot(py,px,pn);
                 xl=1:nc;
                plot(xl,stb00(xl),'.',xl,stb01(xl),'.',xl,stb02(xl),'.',xl,stb10(xl),'.',xl,stb11(xl),'.',xl,stb20(xl),'.')
                set(gca,'xLim',[0 ni+1]);
                set(gca,'yLim',[-1.5 1.5]);
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

%        % Write imgage volume
%        fname = strcat('./fu.bfloat')
%        fid = fopen(fname,'wb');
%        fwrite(fid,M,'float32');
%        status=fclose(fid);
%        % write image header file
%        hname = strcat('./fu.hdr');
%        fid = fopen(hname,'wb');
%        fprintf(fid,'%d %d %d %d # this line is for xds\n',nx,ny,Rep,nsl);
%        fprintf(fid,'x %d \n',nx);
%        fprintf(fid,'y %d \n',ny);
%        fprintf(fid,'z %d \n',nsl);
%        fprintf(fid,'t %d \n',Rep);
%        fprintf(fid,'byte-order %d \n',1);
%        fprintf(fid,'resolution %d %d %d \n',ResX,ResY,ResZ);
%        fprintf(fid,'#origin %d %d %d \n',0,0,0);
%        status=fclose(fid);
%        %

        unix('if [ -a bfloat ] ; then echo "directory bfloat exists" ; else mkdir bfloat ; fi');

%        % flip dimensions to fit SPM view
        fldata=flipdim(flipdim(M,1),3); 

        for zz=1:nsl
        %    zz
            % copy data
            sldata=squeeze(M(:,:,zz,:)*fftscale); 
            % save file
            eval(sprintf('datfile=''./bfloat/fu_%03.0f.bfloat'';',zz-1));
            fid = fopen(datfile,'wb');
            fwrite(fid,sldata,'float32');
            status=fclose(fid);
            % write header
            eval(sprintf('hdrfile=''./bfloat/fu_%03.0f.hdr'';',zz-1));
            fid = fopen(hdrfile,'wb');
            fprintf(fid,'%d %d %d %d # this line is for xds\n',Nx,Ny,Rep,1);
            fprintf(fid,'x %d \n',Nx);
            fprintf(fid,'y %d \n',Ny);
            fprintf(fid,'z %d \n',1);
            fprintf(fid,'t %d \n',Rep);
            fprintf(fid,'byte-order %d \n',1);
            fprintf(fid,'resolution %g %g %g \n',ResX,ResY,ResZ);
            %fprintf(fid,'origin %d %d %d \n',0,0,0);
            status=fclose(fid);
        end
        
       % write b-header
        eval(sprintf('bhdrfile=''./bfloat/fu.bhdr'';'));
        fid = fopen(bhdrfile,'wb');
        fprintf(fid,'          cols: %d \n', nx );
        fprintf(fid,'          rows: %d \n', ny );
        fprintf(fid,'       nslices: %d \n', nsl );
        fprintf(fid,' n_time_points: %d \n', Rep );
        fprintf(fid,'   slice_thick: %d \n', dSlice );
        fprintf(fid,'    top_left_r: %g \n', -FovX/2 );
        fprintf(fid,'    top_left_a: %g \n', -FovY/2 );
        fprintf(fid,'    top_left_s: %g \n', -nsl*dSlice/2 );
        fprintf(fid,'   top_right_r: %g \n', FovX/2 );
        fprintf(fid,'   top_right_a: %g \n', -FovY/2 );
        fprintf(fid,'   top_right_s: %g \n', -nsl*dSlice/2 );
        fprintf(fid,'bottom_right_r: %g \n', FovX/2 );
        fprintf(fid,'bottom_right_a: %g \n', FovY/2 );
        fprintf(fid,'bottom_right_s: %g \n', -nsl*dSlice/2 );
        fprintf(fid,'      normal_r: %d \n', 0 );
        fprintf(fid,'      normal_a: %d \n', 0 );
        fprintf(fid,'      normal_s: %d \n', 1 );
        fprintf(fid,'      image_te: %d \n', TE );
        fprintf(fid,'      image_tr: %d \n', TR );
        fprintf(fid,'      image_ti: %d \n', TE );
        fprintf(fid,'    flip_angle: %g \n',1.5708 );
        status=fclose(fid);
        

    end % dowrite
