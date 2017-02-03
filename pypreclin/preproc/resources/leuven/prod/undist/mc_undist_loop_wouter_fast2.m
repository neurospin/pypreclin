function [W] = mc_undist_loop_wouter_fast2(V)

nx=V.nx;
ny=V.ny;

tim1=V.im1';
tim2=V.im2';

mask_V = tim1;
mask_V(V.tmask' <= 0) = NaN;

indx_mask = ~isnan(mask_V);
tim1_mask = tim1(indx_mask);
    
b00_best=0;
b01_best=0;
b02_best=0;
b10_best=0;
b11_best=0;
b20_best=0;

indx = ~isnan(tim2+mask_V);
tim2r_OK = tim2(indx);
tim1_OK = tim1(indx);
costt = sum((tim2r_OK-tim1_OK).^2);

tim2r_best = tim2;
bim_best = zeros(size(tim2r_best));

if(costt == 0)
	cost_best = 0;
	scale = 0;
	iter = 0;
else
	costb = costt + 2*sum(tim2r_OK.*tim1_OK) +1;
	cost_best = costt/costb;
	
	b01_adj_best = 0;
	b02_adj_best = 0;
	b10_adj_best = 0;
	b11_adj_best = 0;
	b20_adj_best = 0;
	
	X = (1:nx)' * ones(1,ny);
	adjustX_lin = (((1:nx)'-(nx+1)/2)/((nx+1)/2)) * ones(1,ny);
	adjustX_sqr = adjustX_lin .^2;
	adjustY_lin = ones(nx,1) * (((1:ny)-(ny+1)/2)/((ny+1)/2));
	adjustY_sqr = adjustY_lin .^2;
	adjustXYsqr = adjustX_lin.*adjustY_lin;
	
	for scale = V.scale_list
		iter  = 0;
		changed = 1;
		
		while changed
			iter = iter+1;
			changed = 0;
			
			if scale ~= 0
				db00list = V.db.list1 * V.db.step(1) * scale;
				db01list = V.db.list2 * V.db.step(2) * scale;
				db02list = V.db.list3 * V.db.step(3) * scale;
				db10list = V.db.list4 * V.db.step(4) * scale;
				db11list = V.db.list5 * V.db.step(5) * scale;
				db20list = V.db.list6 * V.db.step(6) * scale;
				
				b00 = b00_best;
				b01 = b01_best;
				b02 = b02_best;
				b10 = b10_best;
				b11 = b11_best;
				b20 = b20_best;
				
				b01_adj = b01_adj_best;
				b02_adj = b02_adj_best;
				b10_adj = b10_adj_best;
				b11_adj = b11_adj_best;
				b20_adj = b20_adj_best;
				
% % 				bimx = X + b00 + b01_adj + b02_adj + b10_adj + b11_adj + b20_adj;
				bimx = X + bim_best;
				
				for db00 = db00list
					if db00 == 0
						for db01 = db01list
							if db01 == 0
								for db02 = db02list
									if db02 == 0
										for db10 = db10list
											if db10 == 0
												for db11 = db11list
													if db11 == 0
														for db20 = db20list
															if(db20 ~= 0)
																b20 = b20_best + db20;
																%
																bimx = bimx-b20_adj;
																b20_adj = b20 .* adjustX_sqr ;
																bimx = bimx + b20_adj;
																%
																tim2r = n_interp1_mc(tim2,bimx);
																%
																tim2r_OK = tim2r(indx_mask);
																costt = sum((tim2r_OK-tim1_mask).^2);
																if(isnan(costt))
																	indx = ~isnan(tim2r+mask_V);
																	tim2r_OK = tim2r(indx);
																	tim1_OK = tim1(indx);
																	costt = sum((tim2r_OK-tim1_OK).^2);
																	% % costb = sum(tim2r_OK.^2+tim1_OK.^2) +1;
																	costb = costt + 2*sum(tim2r_OK.*tim1_OK) +1;
																else
																	% % costb = sum(tim2r_OK.^2+tim1_mask.^2) +1;
																	costb = costt + 2*sum(tim2r_OK.*tim1_mask) +1;
																end;
																cost=costt/costb;
																%
																if cost<cost_best
																	cost_best = cost;
																	changed = 1;
																	b00_best = b00;
																	b01_best = b01;
																	b02_best = b02;
																	b10_best = b10;
																	b11_best = b11;
																	b20_best = b20;
																	b01_adj_best = b01_adj;
																	b02_adj_best = b02_adj;
																	b10_adj_best = b10_adj;
																	b11_adj_best = b11_adj;
																	b20_adj_best = b20_adj;
																	tim2r_best = tim2r;
																	bim_best = bimx - X;
																end
															end
														end
													else
														b11 = b11_best + db11;
														
														bimx = bimx - b11_adj;
														b11_adj = b11*adjustXYsqr;
														bimx = bimx + b11_adj;
														
														tim2r = n_interp1_mc(tim2,bimx);
														
														tim2r_OK = tim2r(indx_mask);
														costt = sum((tim2r_OK-tim1_mask).^2);
														if(isnan(costt))
															indx = ~isnan(tim2r+mask_V);
															tim2r_OK = tim2r(indx);
															tim1_OK = tim1(indx);
															costt = sum((tim2r_OK-tim1_OK).^2);
															% % costb = sum(tim2r_OK.^2+tim1_OK.^2) +1;
															costb = costt + 2*sum(tim2r_OK.*tim1_OK) +1;
														else
															% % costb = sum(tim2r_OK.^2+tim1_mask.^2) +1;
															costb = costt + 2*sum(tim2r_OK.*tim1_mask) +1;
														end;
														cost=costt/costb;
														%
														if cost<cost_best
															cost_best = cost;
															changed = 1;
															b00_best = b00;
															b01_best = b01;
															b02_best = b02;
															b10_best = b10;
															b11_best = b11;
															b20_best = b20;
															b01_adj_best = b01_adj;
															b02_adj_best = b02_adj;
															b10_adj_best = b10_adj;
															b11_adj_best = b11_adj;
															b20_adj_best = b20_adj;
															tim2r_best = tim2r;
															bim_best = bimx - X;
														end
													end
												end
											else
												b10 = b10_best + db10;
												
												bimx = bimx - b10_adj;
												b10_adj = b10 * adjustX_lin;
												bimx = bimx + b10_adj;
												
												tim2r = n_interp1_mc(tim2,bimx);
												
												tim2r_OK = tim2r(indx_mask);
												costt = sum((tim2r_OK-tim1_mask).^2);
												if(isnan(costt))
													indx = ~isnan(tim2r+mask_V);
													tim2r_OK = tim2r(indx);
													tim1_OK = tim1(indx);
													costt = sum((tim2r_OK-tim1_OK).^2);
													% % costb = sum(tim2r_OK.^2+tim1_OK.^2) +1;
													costb = costt + 2*sum(tim2r_OK.*tim1_OK) +1;
												else
													% % costb = sum(tim2r_OK.^2+tim1_mask.^2) +1;
													costb = costt + 2*sum(tim2r_OK.*tim1_mask) +1;
												end;
												cost=costt/costb;
												%
												if cost<cost_best
													cost_best = cost;
													changed = 1;
													b00_best = b00;
													b01_best = b01;
													b02_best = b02;
													b10_best = b10;
													b11_best = b11;
													b20_best = b20;
													b01_adj_best = b01_adj;
													b02_adj_best = b02_adj;
													b10_adj_best = b10_adj;
													b11_adj_best = b11_adj;
													b20_adj_best = b20_adj;
													tim2r_best = tim2r;
													bim_best = bimx - X;
												end
											end
										end
									else
										b02 = b02_best + db02;
										
										bimx = bimx - b02_adj;
										b02_adj = b02*adjustY_sqr;
										bimx = bimx + b02_adj;
										
										tim2r = n_interp1_mc(tim2,bimx);
										
										tim2r_OK = tim2r(indx_mask);
										costt = sum((tim2r_OK-tim1_mask).^2);
										if(isnan(costt))
											indx = ~isnan(tim2r+mask_V);
											tim2r_OK = tim2r(indx);
											tim1_OK = tim1(indx);
											costt = sum((tim2r_OK-tim1_OK).^2);
											% % costb = sum(tim2r_OK.^2+tim1_OK.^2) +1;
											costb = costt + 2*sum(tim2r_OK.*tim1_OK) +1;
										else
											% % costb = sum(tim2r_OK.^2+tim1_mask.^2) +1;
											costb = costt + 2*sum(tim2r_OK.*tim1_mask) +1;
										end;
										cost=costt/costb;
										%
										if cost<cost_best
											cost_best = cost;
											changed = 1;
											b00_best = b00;
											b01_best = b01;
											b02_best = b02;
											b10_best = b10;
											b11_best = b11;
											b20_best = b20;
											b01_adj_best = b01_adj;
											b02_adj_best = b02_adj;
											b10_adj_best = b10_adj;
											b11_adj_best = b11_adj;
											b20_adj_best = b20_adj;
											tim2r_best = tim2r;
											bim_best = bimx - X;
										end
									end
								end
							else
								b01 = b01_best + db01;
								
								bimx = bimx - b01_adj;
								b01_adj = b01*adjustY_lin;
								bimx = bimx + b01_adj;
								
								tim2r = n_interp1_mc(tim2,bimx);
								
								tim2r_OK = tim2r(indx_mask);
								costt = sum((tim2r_OK-tim1_mask).^2);
								if(isnan(costt))
									indx = ~isnan(tim2r+mask_V);
									tim2r_OK = tim2r(indx);
									tim1_OK = tim1(indx);
									costt = sum((tim2r_OK-tim1_OK).^2);
									% % costb = sum(tim2r_OK.^2+tim1_OK.^2) +1;
									costb = costt + 2*sum(tim2r_OK.*tim1_OK) +1;
								else
									% % costb = sum(tim2r_OK.^2+tim1_mask.^2) +1;
									costb = costt + 2*sum(tim2r_OK.*tim1_mask) +1;
								end;
								cost=costt/costb;
								%
								if cost<cost_best
									cost_best = cost;
									changed = 1;
									b00_best = b00;
									b01_best = b01;
									b02_best = b02;
									b10_best = b10;
									b11_best = b11;
									b20_best = b20;
									b01_adj_best = b01_adj;
									b02_adj_best = b02_adj;
									b10_adj_best = b10_adj;
									b11_adj_best = b11_adj;
									b20_adj_best = b20_adj;
									tim2r_best = tim2r;
									bim_best = bimx - X;
								end
							end
						end
					else
						bimx = bimx - b00;
						b00 = b00_best + db00;
						bimx = bimx + b00;
						
						tim2r = n_interp1_mc(tim2,bimx);
						
						tim2r_OK = tim2r(indx_mask);
						costt = sum((tim2r_OK-tim1_mask).^2);
						if(isnan(costt))
							indx = ~isnan(tim2r+mask_V);
							tim2r_OK = tim2r(indx);
							tim1_OK = tim1(indx);
							costt = sum((tim2r_OK-tim1_OK).^2);
							% % costb = sum(tim2r_OK.^2+tim1_OK.^2) +1;
							costb = costt + 2*sum(tim2r_OK.*tim1_OK) +1;
						else
							% % costb = sum(tim2r_OK.^2+tim1_mask.^2) +1;
							costb = costt + 2*sum(tim2r_OK.*tim1_mask) +1;
						end;
						cost=costt/costb;
						%
						if cost<cost_best
							cost_best = cost;
							changed = 1;
							b00_best = b00;
							b01_best = b01;
							b02_best = b02;
							b10_best = b10;
							b11_best = b11;
							b20_best = b20;
							b01_adj_best = b01_adj;
							b02_adj_best = b02_adj;
							b10_adj_best = b10_adj;
							b11_adj_best = b11_adj;
							b20_adj_best = b20_adj;
							tim2r_best = tim2r;
							bim_best = bimx - X;
						end
					end
				end
			end
		end
	end
end
tim2r_best(isnan(tim2r_best)) = 0;

% Output structure
W.im3 = tim2r_best';
W.bim = bim_best';
W.b00 = b00_best;
W.b01 = b01_best;
W.b02 = b02_best;
W.b10 = b10_best;
W.b11 = b11_best;
W.b20 = b20_best;
W.cst = cost_best;
W.scale = scale;
W.iter = iter;
