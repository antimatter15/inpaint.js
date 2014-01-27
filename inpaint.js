// This file isn't very well documented, but that shouldn't be a huge problem
// because most of this is a fairly straightforward port of the scikit-image
// implementation which can be found in 
// https://github.com/chintak/scikit-image/blob/inpaint/skimage/filter/_inpaint_fmm.pyx

function InpaintTelea(image, mask, radius){
	if(!radius) radius = 5;

	var LARGE_VALUE = 1e6;
	var SMALL_VALUE = 1e-6;

	var outside = new jsfeat.matrix_t(mask.cols, mask.rows, jsfeat.U8C1_t);
	var flag = new jsfeat.matrix_t(mask.cols, mask.rows, jsfeat.U8C1_t);
	var u = new jsfeat.matrix_t(mask.cols, mask.rows, jsfeat.F32_t | jsfeat.C1_t);

	for(var i = 0; i < mask.data.length; i++){
		if(!mask.data[i]) continue;
		// this is the equivalent of doing a morphological dilation with
		// a 1-pixel cross structuring element
		outside.data[i + 1] = outside.data[i] = outside.data[i - 1] = outside.data[i + outside.cols] = outside.data[i - outside.cols] = 1
	}
	
	for(var i = 0; i < mask.data.length; i++){
		flag.data[i] = (outside.data[i] * 2) - (mask.data[i] ^ outside.data[i])
		if(flag.data[i] == 2) // UNKNOWN
			u.data[i] = LARGE_VALUE;
	}

	function eikonal(n1, n2){
		var u1 = u.data[n1],
			u2 = u.data[n2];
		if(flag.data[n1] == 0 /*KNOWN*/){
			if(flag.data[n2] == 0 /*KNOWN*/){
				var perp = Math.sqrt(2 - Math.pow(u1 - u2, 2)); // perpendicular distance
				var s = (u1 + u2 - perp) * 0.5; // average distance
				if(s >= u1 && s >= u2){
					return s;
				}else{
					s += perp;
					if(s >= u1 && s >= u2){
						return s;
					}
				}
			}else{
				return 1 + u1
			}
		}else if(flag.data[n2] == 0 /*KNOWN*/){
			return 1 + u2;
		}
		return LARGE_VALUE
	}

	function inpaint_point(n){
		var h = flag.rows, w = flag.cols;
		var Ia = 0, Jx = 0, Jy = 0, norm = 0;
		var gradx_u = grad_func(u, n, 1),
			grady_u = grad_func(u, n, flag.cols); 
		
		var i = n % flag.cols,
			j = Math.floor(n / flag.cols);


		for(var k = 0; k < indices_centered.length; k++){
			var nb = n + indices_centered[k];
			var i_nb = nb % flag.cols,
				j_nb = Math.floor(nb / flag.cols);

			if(i_nb <= 1 || j_nb <= 1 || i_nb >= w - 1 || j_nb >= h - 1) continue;

			if(flag.data[nb] != 0 /*KNOWN*/) continue; 

			var rx = i - i_nb,
				ry = j - j_nb;

			var geometric_dst = 1 / ((rx * rx + ry * ry) * Math.sqrt(rx * rx + ry * ry))
			var levelset_dst = 1 / (1 + Math.abs(u.data[nb] - u.data[n]))
			var direction = Math.abs(rx * gradx_u + ry * grady_u) + SMALL_VALUE;
			var weight = geometric_dst * levelset_dst * direction;
			var gradx_img = grad_func(image, nb, 1),
				grady_img = grad_func(image, nb, flag.cols);
			
			Ia += weight * image.data[nb]
			Jx -= weight * gradx_img * rx
			Jy -= weight * grady_img * ry
			norm += weight
			
		}
		image.data[n] = Ia / norm + (Jx + Jy) / Math.sqrt(Jx * Jx + Jy * Jy);			
	}


	function grad_func(array, n, step){
		if(flag.data[n + step] != 2 /* UNKNOWN */){
			if(flag.data[n - step] != 2){
				return array.data[n + step] - array.data[n - step] * 0.5
			}else{
				return array.data[n + step] - array.data[n]
			}
		}else{
			if(flag.data[n - step] != 2){
				return array.data[n] - array.data[n - step]
			}else{
				return 0
			}
		}
	}

	var heap = new HeapQueue(function(a, b){
		return a[0] - b[0] // sort by first thingy
	})
	
	for(var i = 0; i < mask.data.length; i++){
		if(flag.data[i] == 1) // BAND
			heap.push([u.data[i], i]);
	}
	
	var indices_centered = []
	// generate a mask for a circular structuring element
	for(var i = -radius; i <= radius; i++){
		var h = Math.floor(Math.sqrt(radius * radius - i * i))
		for(var j = -h; j <= h; j++)
			indices_centered.push(i + j * flag.cols);
	}

	while(heap.length){
		var n = heap.pop()[1];
		var i = n % mask.cols,
			j = Math.floor(n / mask.cols);
		
		flag.data[n] = 0; // KNOWN

		if(i <= 1 || j <= 1 || i >= mask.cols - 1 || j >= mask.rows - 1) continue;

		for(var k = 0; k < 4; k++){
			var nb = n + [-mask.cols, -1, mask.cols, 1][k];
			if(flag.data[nb] != 0){ // not KNOWN

				u.data[nb] = Math.min(eikonal(nb - mask.cols, nb - 1),
                                      eikonal(nb + mask.cols, nb - 1),
                                      eikonal(nb - mask.cols, nb + 1),
                                      eikonal(nb + mask.cols, nb + 1))

				if(flag.data[nb] == 2){ // UNKNOWN
					flag.data[nb] = 1; // BAND
					heap.push([u.data[nb], nb])
					inpaint_point(nb)
				}
			}
		}
	}
}



