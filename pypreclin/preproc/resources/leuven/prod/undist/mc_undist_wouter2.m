function [] = mc_undist_wouter2(edir)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home;
fprintf('Motion Compensation Undistortion: ''mc_undist_wouter2.m''\n');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seq=100101;
% % seq=110111;
% % seq=111111;
options=0;
PE=0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg = sprintf('------ Undistortion Parameters ------\n'); disp(msg);
%msg = sprintf(' Target Dir          = %d\n',tdir); disp(msg);
msg = sprintf(' UD sequence         = %d\n',seq); disp(msg);
msg = sprintf(' UD options          = %d\n',options); disp(msg);
%msg = sprintf(' FFT Scale Factor    = %g\n',fftscale); disp(msg);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dd=dir([edir '*_msk*.nii']);
path_mask=[edir dd(1).name];

V_mask=spm_vol(path_mask);

mask = spm_read_vols(V_mask);

dd2=dir([edir '*_epi*.nii']);
path_epi=[edir dd2(1).name];

V_epi=spm_vol(path_epi);
xdim = V_epi(1).dim(1);
ydim = V_epi(1).dim(2);
zdim = V_epi(1).dim(3);
tdim = length(V_epi);

switch(spm_type(V_epi(1).dt(1)))
	case 'float64'
		dtype = 'double';
		% % mask_smooth = smooth3(double(mask > 0),'box',7);
		mask_smooth = sqrt(smooth3(double(smooth3(double(mask > 0), 'gaussian', 7, 1) > 0.25),'gaussian', 7, 1));
 	case 'float32'
		dtype = 'single';
		% % mask_smooth = smooth3(single(mask > 0),'box',7);
		mask_smooth = sqrt(smooth3(single(smooth3(single(mask > 0), 'gaussian', 7, 1) > 0.25),'gaussian', 7, 1));
	otherwise
		dtype = 'single';
		% % mask_smooth = smooth3(single(mask > 0),'box',7);
		mask_smooth = sqrt(smooth3(single(smooth3(single(mask > 0), 'gaussian', 7, 1) > 0.25),'gaussian', 7, 1));
end;

dx = round(xdim/10); dy = round(ydim/10); dz = round(zdim/10);
select_corners = zeros([xdim ydim zdim]);
select_corners(1:dx, 1:dy, 1:dz) = 1;
select_corners(end-dx:end, 1:dy, 1:dz) = 1;
select_corners(1:dx, end-dx:end, 1:dz) = 1;
select_corners(end-dx:end, end-dx:end, 1:dz) = 1;
select_corners(1:dx, 1:dy, end-dz:end) = 1;
select_corners(end-dx:end, 1:dy, end-dz:end) = 1;
select_corners(1:dx, end-dy:end, end-dz:end) = 1;
select_corners(end-dx:end, end-dy:end, end-dz:end) = 1;
select_corners = select_corners > 0;

% Reading piece by piece, allows for progress 'bar' (Slow spm_read_vols on
% Mac OS X)
% T: Volumes 'preprocessed' (masked) for alignment
T = zeros([xdim ydim zdim tdim], dtype);
% M: Raw volumes / output
M = zeros([xdim ydim zdim tdim], dtype);
% Mean volume
Tm = zeros([xdim ydim zdim], dtype);

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fprintf('\n');
% % fprintf(repmat(' ', [1 80]));
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:tdim
	
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 	fprintf(repmat('\b', [1 80]));
% % 	infotxt = sprintf('reading volume %5s of %5s', ...
% % 		sprintf('%d', n), sprintf('%d', tdim));
% % 	fprintf('%80s', infotxt);
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	vol = spm_read_vols(V_epi(n));
	M(:,:,:,n) = vol; % store raw volume
	Tm = Tm + vol;
	
	noise = max(vol(select_corners));
	vol = vol.*mask_smooth;
	vol(vol<noise) = 0;
	T(:,:,:,n) = vol; % 'prepared' volume vor alignment
end;
Tm = Tm / tdim;

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fprintf(repmat('\b', [1 80]));
% % fprintf('%80s\n', 'reading volumes ... done');
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine target volume
[tmp vox] = max(Tm(:).*mask_smooth(:));
[xv yv zv] = ind2sub([xdim ydim zdim], vox);
sig = smooth(detrend(permute(M(xv, yv, zv,:), [4 1 2 3])).^2);
thres = prctile(sig,10);
vect=find(sig<thres);

res = ones(1,length(vect));
for n=1:length(vect)
	Tr=abs(M(:,:,:,vect(n)).*mask_smooth-Tm);
	res(n)=mean(Tr(:));
end
[tmp nt] = min(res);
nt = vect(nt);
Tt = T(:,:,:,nt);
% % Tt = Tm.*mask_smooth;

%
if options==0
	dowrite=1 ;
	doquad=1;
elseif options==1
	dowrite=1 ;
	doquad=1;
elseif options==2 % same as options=0 but 4 quadrants
	dowrite=1 ;
	doquad=4;
elseif options==3 % same as options=1 but 4 quadrants
	dowrite=1 ;
	doquad=4;
end

if PE==2||PE==0
	% temp output array
	rnx = xdim;
	rny = ydim;
elseif PE==1||PE==-1
	rnx = ydim;
	rny = xdim;
end
% % [X,Y] = meshgrid(1:rny,1:rnx);
X = ones([rnx, 1], dtype) * (1:rny);
Y = (1:rnx)' * ones([1, rny], dtype);

% smoothing
% % psf=fspecial('gaussian',5,1);
psf=fspecial('gaussian',5,sqrt(2)/2);
% % psf=fspecial('gaussian',7,sqrt(2));

% init statistics vars

% %%%%%%%%
% slice loop start

for slc=1:zdim
% % for slc=1:10
	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	msg = sprintf('undist. slice : %d',slc); disp(msg);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Target Image
	im1 = Tt(:,:,slc);
% % 	im1 = Tm(:,:,slc);
	mask_slice = mask_smooth(:,:,slc);
	switch(PE)
		% case 1
		%	im1 = im1;
		case -1
			im1 = rot90(im1, 1);
			mask_slice = rot90(mask_slice, 1);
		case 0
			im1 = rot90(im1, 2);
			mask_slice = rot90(mask_slice, 2);
		case 1
			im1 = rot90(im1, 3);
			mask_slice = rot90(mask_slice, 3);
	end;
	
	% %         % smoothed target image
	% %         im=imfilter(im1,psf,'same','conv');
	% %         tim1=deconvlucy(edgetaper(im,psf),psf,50);
	tim1=imfilter(im1,psf,'same','conv', 'replicate');
		
	%
	qnn=doquad;
	stb00=zeros(tdim,qnn);
	stb01=zeros(tdim,qnn);
	stb02=zeros(tdim,qnn);
	stb10=zeros(tdim,qnn);
	stb11=zeros(tdim,qnn);
	stb20=zeros(tdim,qnn);
	
	% init structure for distortion correction
	V.tmask = double(mask_slice);
	V.nx=rnx;
	V.ny=rny;
	
	% %%%%%%%%%
	% time loop start
	
	for nc=1:tdim
		
		% align data images
		im2 = T(:,:,slc,nc);
		im2_raw = M(:,:,slc,nc);
		switch(PE)
			% case 2
			%	im2 = im2;
			case -1
				im2 = rot90(im2, 1);
				im2_raw = rot90(im2_raw, 1);
			case 0
				im2 = rot90(im2, 2);
				im2_raw = rot90(im2_raw, 2);
			case 1
				im2 = rot90(im2, 3);
				im2_raw = rot90(im2_raw, 3);
		end;
		% smoothed distorted image
		% %             im=imfilter(im2,psf,'same','conv');
		% %             tim2=deconvlucy(edgetaper(im,psf),psf,20);
		tim2=imfilter(im2,psf,'same','conv', 'replicate');
		
		% init final image
		im3=zeros(size(im2));
		
		W.bim=zeros(size(im1));
		
		qim1 = tim1;
		qim2 = tim2;
		
		Wqim=zeros(size(im1));
		
		if doquad==4
			% qim1=zeros(size(tim1,1),size(tim1,2),4);
			qim1 = zeros([size(tim1) 4]);
			qim1(1:end/2,1:end/2,1) = tim1(1:end/2,1:end/2);
			qim1(1:end/2,end/2+1:end,2) = tim1(1:end/2,end/2+1:end);
			qim1(end/2+1:end,end/2+1:end,3) = tim1(end/2+1:end,end/2+1:end);
			qim1(end/2+1:end,1:end/2,4) = tim1(end/2+1:end,1:end/2);
			% qim2=zeros(size(tim2,1),size(tim2,2),4);
			qim2 = zeros([size(tim2) 4]);
			qim2(1:end/2,1:end/2,1) = tim2(1:end/2,1:end/2);
			qim2(1:end/2,end/2+1:end,2) = tim2(1:end/2,end/2+1:end);
			qim2(end/2+1:end,end/2+1:end,3) = tim2(end/2+1:end,end/2+1:end);
			qim2(end/2+1:end,1:end/2,4) = tim2(end/2+1:end,1:end/2);
		end
		
		for qq=1:qnn
			
			% % 				if (length(find(tmask))>100);
			% % 				if (sum(im1(:)+im2(:) > 0) > 0.01*xdim*ydim)
			if (((sum(im1(:))/sum(im2(:)) > 0.8) && (sum(im1(:))/sum(im2(:)) < 1.25)) && sum(im1(:)+im2(:) > 0) > 100)
				
				tqim1=qim1(:,:,qq);
				tqim2=qim2(:,:,qq);
				
				% target image
				V.im1 = double(tqim1 .* tim1);
				
				
				% distorted image
				V.im2 = double(tqim2 .* tim2);
				
				% %% vars for optimization loop				
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
				
				V.db.list1 = 0; if V.db.dolist(1); V.db.list1=[-1 0 1]; end
				V.db.list2 = 0; if V.db.dolist(2); V.db.list2=[-1 0 1]; end
				V.db.list3 = 0; if V.db.dolist(3); V.db.list3=[-1 0 1]; end
				V.db.list4 = 0; if V.db.dolist(4); V.db.list4=[-1 0 1]; end
				V.db.list5 = 0; if V.db.dolist(5); V.db.list5=[-1 0 1]; end
				V.db.list6 = 0; if V.db.dolist(6); V.db.list6=[-1 0 1]; end
				
				
				% %                  % subroutine
				W=mc_undist_loop_wouter(V) ; % updated (28/06/2011)
				
				% ouput
				stb00(nc,qq)=W.b00;
				stb01(nc,qq)=W.b01;
				stb02(nc,qq)=W.b02;
				stb10(nc,qq)=W.b10;
				stb11(nc,qq)=W.b11;
				stb20(nc,qq)=W.b20;
				scale=W.scale;
				iter=W.iter;
				
				if options==1;
					fprintf(1,'rep = %d: scale = %d (%d): b00 = %1.3f, b01 = %1.3f, b02 = %1.3f, b10 = %1.3f, b11 = %1.3f, b20 = %1.3f, cost = %1.5f\n',nc,scale,iter,W.b00,W.b01,W.b02,W.b10,W.b11,W.b20,W.cst*1e16);
				end
				
				if(any(W.bim(:) ~= 0))
					qim3 = interp2(im2_raw,X+W.bim,Y,'cubic');
% % 					qim3 = n_interp1_mc(im2_raw',(X+W.bim)')';
				else
					qim3 = im2_raw;
				end;
				
			else % tmask
				
% % 				qim3=im1;
				qim3=im2_raw;
				
				if options==1;
					fprintf(1,'rep = %d\n',nc);
				end
			end % tmask
			
			if doquad==1
				im3=qim3;
			elseif doquad==4
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
		
		switch(PE)
			case 2
				M(:,:,slc,nc) = im3;
			case -1
				M(:,:,slc,nc) = rot90(im3, 3);
			case 0
				M(:,:,slc,nc) = rot90(im3, 2);
			case +1
				M(:,:,slc,nc) = rot90(im3, 1);
		end;

	end % nc
	
end % slc

if dowrite==1
	fname=[edir 'u_' dd2(1).name];
	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf('Writing results to disk: ''%s''.\n', fname);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	for n=1:length(V_epi)
		V_out = V_epi(n);
		V_out.fname = fname;
		spm_write_vol(V_out, M(:,:,:,n).*mask_smooth);
	end
end % dowrite

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Motion Compensation Undistortion: ... done\n');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
