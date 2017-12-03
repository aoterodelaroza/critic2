    // grid and Cartesian axes
    // int nmaxgrid = (int) (ceil(c2::scenerad * AUTOANG)+1e-5);
    // vec3 r1 = {}, r2 = {};
    // vec4 rgb = {1.f,1.f,1.f,1.f};
    // rgb.x = 1.f; rgb.y = 0.f; rgb.z = 0.f;
    // r1.y = 0.f; r1.x = -nmaxgrid/AUTOANG;
    // r2.y = 0.f; r2.x =  nmaxgrid/AUTOANG;
    // drawCylinder(r1,r2,radgrid,rgb,0,false);
    // rgb.x = 0.f; rgb.y = 1.f; rgb.z = 0.f;
    // r1.x = 0.f; r1.y = -nmaxgrid/AUTOANG;
    // r2.x = 0.f; r2.y =  nmaxgrid/AUTOANG;
    // drawCylinder(r1,r2,radgrid,rgb,0,false);
    // rgb.x = 0.f; rgb.y = 0.f; rgb.z = 1.f;
    // r1.y = 0.f; r1.z = -nmaxgrid/AUTOANG;
    // r2.y = 0.f; r2.z =  nmaxgrid/AUTOANG;
    // drawCylinder(r1,r2,radgrid,rgb,0,false);
    // for (int i = 1; i <= nmaxgrid; i++){
    //   r1.z = 0.f; r2.z = 0.f;
    //   rgb.x = 1.f; rgb.y = 1.f; rgb.z = 1.f; rgb.w = 1.f; 
    //   r1.x = i/AUTOANG; r1.y = -nmaxgrid/AUTOANG;
    //   r2.x = i/AUTOANG; r2.y =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    //   r1.x = -i/AUTOANG; r1.y = -nmaxgrid/AUTOANG;
    //   r2.x = -i/AUTOANG; r2.y =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    //   r1.y =  i/AUTOANG; r1.x = -nmaxgrid/AUTOANG;
    //   r2.y =  i/AUTOANG; r2.x =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    //   r1.y = -i/AUTOANG; r1.x = -nmaxgrid/AUTOANG;
    //   r2.y = -i/AUTOANG; r2.x =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    // }

    // rotation sphere
    // if (llock){
    //   vec3 v0;
    //   float rad0 = 0.25f * v_pos[2];
    //   float rad1;
    //   vec4 rgb;
    //   vec3 uvec = vec3(0.f,0.f,1.f);

    //   v0 = inverse(mat3(m_world)) * uvec;
    //   v0 = v0 * rad0;
    //   rgb = {0.6f,0.2f,0.8f,0.8f};
    //   rad1 = 0.05f * rad0;
    //   drawSphere(v0,rad1,rgb,3,false);

    //   v0 = inverse(mat3(crot0_l)) * uvec;
    //   v0 = v0 * rad0;
    //   rgb = {1.0f,0.8f,0.0f,0.8f};
    //   rad1 = 0.05f * rad0;
    //   drawSphere(v0,rad1,rgb,3,false);

    //   v0 = vec3(0.f,0.f,0.f);
    //   rgb = {1.0f,1.0f,1.0f,0.4f};
    //   drawSphere(v0,rad0,rgb,3,true);
    // }
