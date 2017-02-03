function [W] = mc_undist_loop_test(V)

nx=V.nx;
ny=V.ny;

tim1=V.im1';
tim2=V.im2';
    
b00_best=0;
b01_best=0;
b02_best=0;
b10_best=0;
b11_best=0;
b20_best=0;

adj01_best = 0;
adj02_best = 0;
adj10_best = 0;
adj11_best = 0;
adj20_best = 0;

mask_V = V.tmask' > 0;
tim1_mask = tim1(mask_V);
[X,Y] = meshgrid(1:nx,1:ny);
X = X';

tim2r_OK = tim2(mask_V);
costt = sum((tim2r_OK-tim1_mask).^2);
if(~isnan(costt))
	costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
else
	indx = mask_V & ~isnan(tim2+tim1);
	tim2r_OK = tim2(indx);
	tim1_OK = tim1(indx);
	costt = sum((tim2r_OK-tim1_OK).^2);
	costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
end;
cost_best = costt / costb;
tim2r_best = tim2;
bim_best = 0;
scale = 0;
iter = 0;
% % cost_best = 9e99;

adjustX_lin = (((1:nx)'-(nx+1)/2)/((nx+1)/2)) * ones(1,ny);
adjustX_sqr = adjustX_lin .^2;
adjustY_lin = ones(nx,1)*(((1:ny)-(ny+1)/2)/((ny+1)/2));
adjustY_sqr = adjustY_lin .^2;
adjustXYsqr = adjustY_lin .* adjustX_lin;

for scale = V.scale_list
	iter  = 0;
	changed = 2;
	
	while (changed > 0 && iter < 25)
		iter = iter+1;
		changed = changed - 1;
		
		if scale==0
			db00list = 0;
			db01list = 0;
			db02list = 0;
			db10list = 0;
			db11list = 0;
			db20list = 0;
		else
			db00list = V.db.list1 * V.db.step(1) * scale;
			db01list = V.db.list2 * V.db.step(2) * scale;
			db02list = V.db.list3 * V.db.step(3) * scale;
			db10list = V.db.list4 * V.db.step(4) * scale;
			db11list = V.db.list5 * V.db.step(5) * scale;
			db20list = V.db.list6 * V.db.step(6) * scale;
		end
		
		bim_best = b00_best + adj01_best + adj02_best + adj10_best + adj11_best + adj20_best;
		
		bim_011111_best = bim_best - b00_best;
		bim_101111_best = bim_best - adj01_best;
		bim_110111_best = bim_best - adj02_best;
		bim_111011_best = bim_best - adj10_best;
		bim_111101_best = bim_best - adj11_best;
		bim_111110_best = bim_best - adj20_best;
		
		cost = zeros(numel(db00list),1) + inf;
		bim = cell(numel(db00list),1);
		for n=1:numel(db00list)
			db00 = db00list(n);
			if(db00 ~= 0)
				b00 = b00_best + db00;
				bim = bim_011111_best + b00;
				%
				tim2r = n_interp1_mc(tim2,X+bim);
				%
				tim2r_OK = tim2r(mask_V);
				costt = sum((tim2r_OK-tim1_mask).^2);
				if(~isnan(costt))
					costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
				else
					indx = mask_V & ~isnan(tim2r+tim1);
					tim2r_OK = tim2r(indx);
					tim1_OK = tim1(indx);
					costt = sum((tim2r_OK-tim1_OK).^2);
					costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
				end;
				cost = costt / costb;
			end;
		end;
		%
		if cost<cost_best
			cost_best = cost;
			changed = 2;
			%
			tim2r_best = tim2r;
			b00_best = b00;
			bim_best = bim;
			%
			% % bim_011111_best = bim_best - b00_best;
			bim_101111_best = bim_best - adj01_best;
			bim_110111_best = bim_best - adj02_best;
			bim_111011_best = bim_best - adj10_best;
			bim_111101_best = bim_best - adj11_best;
			bim_111110_best = bim_best - adj20_best;
		end
		
		for db00 = db00list
			if(db00 == 0)
				for db01 = db01list
					if(db01 == 0)
						for db02 = db02list
							if(db02 == 0)
								for db10 = db10list
									if(db10 == 0)
										for db11 = db11list
											if(db11 == 0)
												for db20 = db20list
													if(db20 ~= 0)
													b20 = b20_best + db20;
													%
													bim = bim_111110_best + b20*adjustX_sqr;
													%
													tim2r = n_interp1_mc(tim2,X+bim);
													%
													tim2r_OK = tim2r(mask_V);
													costt = sum((tim2r_OK-tim1_mask).^2);
													if(~isnan(costt))
														costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
													else
														indx = mask_V & ~isnan(tim2r+tim1);
														tim2r_OK = tim2r(indx);
														tim1_OK = tim1(indx);
														costt = sum((tim2r_OK-tim1_OK).^2);
														costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
													end;
													cost = costt / costb;
													%
													if cost<cost_best
														cost_best = cost;
														changed = 1;
														%
														tim2r_best = tim2r;
														%
														bim_best = bim;
														b20_best = b20;
														%
														adj20_best = b20*adjustX_sqr;
														%
														bim_011111_best = bim_best - b00_best;
														bim_101111_best = bim_best - adj01_best;
														bim_110111_best = bim_best - adj02_best;
														bim_111011_best = bim_best - adj10_best;
														bim_111101_best = bim_best - adj11_best;
													end
													end
												end
											else
												b11 = b11_best + db11;
												bim = bim_111101_best + b11*adjustXYsqr;
												%
												tim2r = n_interp1_mc(tim2,X+bim);
												%
												tim2r_OK = tim2r(mask_V);
												costt = sum((tim2r_OK-tim1_mask).^2);
												if(~isnan(costt))
													costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
												else
													indx = mask_V & ~isnan(tim2r+tim1);
													tim2r_OK = tim2r(indx);
													tim1_OK = tim1(indx);
													costt = sum((tim2r_OK-tim1_OK).^2);
													costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
												end;
												cost = costt / costb;
												%
												if cost<cost_best
													cost_best = cost;
													changed = 1;
													%
													tim2r_best = tim2r;
													%
													b11_best = b11;
													adj11_best = b11*adjustXYsqr;
													%
													bim_best = bim;
													%
													bim_011111_best = bim_best - b00_best;
													bim_101111_best = bim_best - adj01_best;
													bim_110111_best = bim_best - adj02_best;
													bim_111011_best = bim_best - adj10_best;
													% % bim_111101_best = bim_best - adj11_best;
													bim_111110_best = bim_best - adj20_best;
												end
											end
										end
									else
										b10 = b10_best + db10;
										bim = bim_111011_best + b10*adjustX_lin;
										%
										tim2r = n_interp1_mc(tim2,X+bim);
										%
										tim2r_OK = tim2r(mask_V);
										costt = sum((tim2r_OK-tim1_mask).^2);
										if(~isnan(costt))
											costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
										else
											indx = mask_V & ~isnan(tim2r+tim1);
											tim2r_OK = tim2r(indx);
											tim1_OK = tim1(indx);
											costt = sum((tim2r_OK-tim1_OK).^2);
											costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
										end;
										cost = costt / costb;
										%
										if cost<cost_best
											cost_best = cost;
											changed = 1;
											%
											tim2r_best = tim2r;
											%
											b10_best = b10;
											adj10_best = b10*adjustX_lin;
											%
											bim_best = bim;
											%
											bim_011111_best = bim_best - b00_best;
											bim_101111_best = bim_best - adj01_best;
											bim_110111_best = bim_best - adj02_best;
											% % bim_111011_best = bim_best - adj10_best;
											bim_111101_best = bim_best - adj11_best;
											bim_111110_best = bim_best - adj20_best;
										end
									end
								end
							else % if(db02 == 0)
								b02 = b02_best + db02;
								bim = bim_110111_best + b02*adjustY_sqr;
								%
								tim2r = n_interp1_mc(tim2,X+bim);
								%
								tim2r_OK = tim2r(mask_V);
								costt = sum((tim2r_OK-tim1_mask).^2);
								if(~isnan(costt))
									costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
								else
									indx = mask_V & ~isnan(tim2r+tim1);
									tim2r_OK = tim2r(indx);
									tim1_OK = tim1(indx);
									costt = sum((tim2r_OK-tim1_OK).^2);
									costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
								end;
								cost = costt / costb;
								%
								if cost<cost_best
									cost_best = cost;
									changed = 1;
									%
									tim2r_best = tim2r;
									%
									b02_best = b02;
									adj02_best = b02*adjustY_sqr;
									%
									bim_best = bim;
									%
									bim_011111_best = bim_best - b00_best;
									bim_101111_best = bim_best - adj01_best;
									% % bim_110111_best = bim_best - adj02_best;
									bim_111011_best = bim_best - adj10_best;
									bim_111101_best = bim_best - adj11_best;
									bim_111110_best = bim_best - adj20_best;
								end
							end % if(db02 == 0)
						end
					else % if(db01 == 0)
						b01 = b01_best + db01;
						bim = bim_101111_best + b01*adjustY_lin;
						%
						tim2r = n_interp1_mc(tim2,X+bim);
						%
						tim2r_OK = tim2r(mask_V);
						costt = sum((tim2r_OK-tim1_mask).^2);
						if(~isnan(costt))
							costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
						else
							indx = mask_V & ~isnan(tim2r+tim1);
							tim2r_OK = tim2r(indx);
							tim1_OK = tim1(indx);
							costt = sum((tim2r_OK-tim1_OK).^2);
							costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
						end;
						cost = costt / costb;
						%
						if cost<cost_best
							cost_best = cost;
							changed = 1;
							%
							tim2r_best = tim2r;
							%
							b01_best = b01;
							adj01_best = b01*adjustY_lin;
							%
							bim_best = bim;
							%
							bim_011111_best = bim_best - b00_best;
							% % bim_101111_best = bim_best - adj01_best;
							bim_110111_best = bim_best - adj02_best;
							bim_111011_best = bim_best - adj10_best;
							bim_111101_best = bim_best - adj11_best;
							bim_111110_best = bim_best - adj20_best;
						end
					end % if(db01 == 0)
				end
			else % if(db00 == 0)
				b00 = b00_best + db00;
				bim = bim_011111_best + b00;
				%
				tim2r = n_interp1_mc(tim2,X+bim);
				%
				tim2r_OK = tim2r(mask_V);
				costt = sum((tim2r_OK-tim1_mask).^2);
				if(~isnan(costt))
					costb = 0.5*sum(tim2r_OK.^2+tim1_mask.^2) + 1;
				else
					indx = mask_V & ~isnan(tim2r+tim1);
					tim2r_OK = tim2r(indx);
					tim1_OK = tim1(indx);
					costt = sum((tim2r_OK-tim1_OK).^2);
					costb = 0.5*sum(tim2r_OK.^2+tim1_OK.^2) + 1;
				end;
				cost = costt / costb;
				%
				if cost<cost_best
					cost_best = cost;
					changed = 1;
					%
					tim2r_best = tim2r;
					%
					b00_best = b00;
					%
					bim_best = bim;
					%
					% % bim_011111_best = bim_best - b00_best;
					bim_101111_best = bim_best - adj01_best;
					bim_110111_best = bim_best - adj02_best;
					bim_111011_best = bim_best - adj10_best;
					bim_111101_best = bim_best - adj11_best;
					bim_111110_best = bim_best - adj20_best;
				end
			end % if(db00 == 0)
		end
	end
end

% % indnan = find(isnan(tim2r_best));
tim2r_best(isnan(tim2r_best)) = 0;

% Output structure
W.im3 = tim2r_best';
W.bim=bim_best';
W.b00 = b00_best;
W.b01 = b01_best;
W.b02 = b02_best;
W.b10 = b10_best;
W.b11 = b11_best;
W.b20 = b20_best;
W.cst = cost_best;
W.scale = scale;
W.iter = iter;
