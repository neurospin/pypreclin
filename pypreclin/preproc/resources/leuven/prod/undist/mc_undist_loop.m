function [W] = mc_undist_loop(V)

nx=V.nx;
ny=V.ny;

tim1=V.im1;
tim2=V.im2;
    
b00_best=0;
b01_best=0;
b02_best=0;
b10_best=0;
b11_best=0;
b20_best=0;

[X,Y] = meshgrid(1:nx,1:ny);

cost_best = 9e99;

adjustX_lin = ones(ny,1) * (([1:nx]-(nx+1)/2)/((nx+1)/2));
adjustX_sqr = adjustX_lin .^2;
adjustY_lin = (([1:ny]'-(ny+1)/2)/((ny+1)/2))*ones(1,nx);
adjustY_sqr = adjustY_lin .^2;
adjust_zero = zeros(ny,nx);

for scale = V.scale_list
    iter  = 0;
    changed = 1;

    while changed
        iter = iter+1;
        changed = 0;

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
 
        for db00 = db00list
            for db01 = db01list
                for db02 = db02list
                    for db10 = db10list
                        for db11 = db11list
                            for db20 = db20list
                                if (db00~=0)+(db01~=0)+(db02~=0)+(db10~=0)+(db11~=0)+(db20~=0)<=1;
                                b00 = b00_best + db00;
                                b01 = b01_best + db01;
                                b02 = b02_best + db02;
                                b10 = b10_best + db10;
                                b11 = b11_best + db11;
                                b20 = b20_best + db20;
                                %                                
                                bim1 = (b00 + b01*adjustY_lin + b02*adjustY_sqr) ;
                                bim2 = (b10 + b11*adjustY_lin) .* adjustX_lin ;
                                bim3 = (b20) .* adjustX_sqr ;
                                bim = bim1 + bim2 + bim3 ;
                                %
                                %tim2r = interp2(tim2,X+bim,Y,'linear');
                                tim2r = bilin_interp(tim2,X+bim,Y);
                                %
                                indx = find( V.tmask & ~(isnan(tim2r)|isnan(tim1)) );
                                tim2r_OK = tim2r(indx);
                                tim1_OK = tim1(indx);
                                costt = (sum(abs(tim2r_OK-tim1_OK).^2));
                                costb = ((sum(tim2r_OK.^2+tim1_OK.^2)/2)+1);
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
                                    tim2r_best = tim2r;
                                    bim_best = bim;
                                end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

indnan = find(isnan(tim2r_best));
tim2r_best(indnan) = 0;

% Output structure
W.im3 = tim2r_best;
W.bim=bim_best;
W.b00 = b00_best;
W.b01 = b01_best;
W.b02 = b02_best;
W.b10 = b10_best;
W.b11 = b11_best;
W.b20 = b20_best;
W.cst = cost_best;
W.scale = scale;
W.iter = iter;
