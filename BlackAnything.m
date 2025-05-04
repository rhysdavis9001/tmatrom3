surf(x(:,:,round(75/2)),y(:,:,round(75/2)),real(totalfield(:,:,round(75/2))))
view([0 90]);
shading interp;
figure(2)
colorbar
title('Total field ')