//JavaScript code for simulation of neutron Laue diffraction pattern at HRC

// 2020/7/28, redefined phiv and phih using Math.atan2
// 2020/7/9, introduced x,y, and z-axes rotation
// 2020/7/7, resolved uy2uz2 division by zero
// 2020/6/25, corrected ux[i]s and expression of reciprocal lattice vectors
// 2020/6/24,  defined 3 reciprocal lattice vectors a*, b* and c* for a sample orientation without rotation (Psi=0) 
// 2020/6/18-19,  introduced lattice constants and sample orientation 
// 2020/6/5
var version = "0.4.1";

var TOFconst = 2.286;       // TOF at 1 m is 2.286/sqrt(E)

var X0 = 0;
var Y0 = 250;
var length1=200;    //unused variable 

var radius=5;

var HD = 40;    // height of center of PSD from incident beam (mm)
var LD = 2800;  // length of PSD (mm)
var L20 = 4000; // distance from sample to PSD in horizontal plane (mm)

var u = new Array(3); // indices, pallarel to the incident beam
var v = new Array(3); // indices, another direction in the horizontal plane including the incidnet beam

var ux = new Array(3);
var vx = new Array(3);

var Rot0 = new Array(3);
var Rot1 = new Array(3);
var Rot2 = new Array(3);
var Rot =[Rot0, Rot1, Rot2];    // 3x3 rotation matrix

// unit vector of primitive translation vectors
var a_unit = new Array(3);
var b_unit = new Array(3);
var c_unit = new Array(3);

// reciprocal lattice vectors
var a_star = new Array(3);
var b_star = new Array(3);
var c_star = new Array(3);

var Hmax;
var Kmax;
var Lmax;

var lambda;             // wavelength for Q-vector
var maxphih = 60.0;     // maximum of the phi angle on the horizontal plane
//var maxphiv = 60.0;
var scaleX=800;
var scaleY=500;

//variables for 3D orientation viewer
var arrow_scale=150;        //arrows for a*, b* and c*: convert A-1 to pixel.
var arrow_HeadLen=20;       //lengths of arrowheads (pixel)
var arrow_HeadWidth=10;     //widths of arrowheads (pixel)
var DetBankAngles=[12.25/180.0*Math.PI, 32.75/180.0*Math.PI, 53.2/180.0*Math.PI, -21.8/180.0*Math.PI];   //angles of the centers of the detector banks (rad)
var DetBankWidth=1300;  // width of the detector banks (mm)
var DetBankScale=0.1;   // convert mm to pixel.


function draw() {
    document.getElementById("verNum").innerHTML=version;
    document.getElementById("verNum2").innerHTML=version;

    set_Lattice();
    draw_DetMap();
    draw_OriViewer();

}

function rot_and_draw(rot_ax_dir) {
    rot_Lattice(rot_ax_dir);
    draw_DetMap();
    draw_OriViewer();
}

function set_Lattice(){

    //input parameters: lattice constants and sample orientation)
    a = Number(document.getElementById('a').value);
    b = Number(document.getElementById('b').value);
    c = Number(document.getElementById('c').value);
    alpha = Number(document.getElementById('alpha').value)/180.0*Math.PI;   // in radian
    beta  = Number(document.getElementById('beta').value)/180.0*Math.PI;    // in radian
    gamma = Number(document.getElementById('gamma').value)/180.0*Math.PI;   // in radian
    u[0] = Number(document.getElementById('u1').value);
    u[1] = Number(document.getElementById('u2').value);
    u[2] = Number(document.getElementById('u3').value);
    v[0] = Number(document.getElementById('v1').value);
    v[1] = Number(document.getElementById('v2').value);
    v[2] = Number(document.getElementById('v3').value);

    // calculation
    let DD = (Math.cos(alpha)-Math.cos(gamma)*Math.cos(beta))/Math.sin(gamma);
    let PP = Math.sqrt(Math.sin(beta)-DD**2.0);

    ux[0] = 2.0*Math.PI*u[0]/a;
    ux[1] = 2.0*Math.PI*(-u[0]/a/Math.tan(gamma)+u[1]/b/Math.sin(gamma));
    ux[2] = 2.0*Math.PI*(u[0]/a*(DD/Math.tan(gamma)-Math.cos(beta))-u[1]/b*DD/Math.sin(gamma)+u[2]/c)/PP;
    vx[0] = 2.0*Math.PI*v[0]/a;
    vx[1] = 2.0*Math.PI*(-v[0]/a/Math.tan(gamma)+v[1]/b/Math.sin(gamma));
    vx[2] = 2.0*Math.PI*(v[0]/a*(DD/Math.tan(gamma)-Math.cos(beta))-v[1]/b*DD/Math.sin(gamma)+v[2]/c)/PP;

    let uy2uz2 = ux[1]**2.0+ux[2]**2.0;    
    let Uabs = Math.sqrt(ux[0]**2.0+uy2uz2);
    let Rvy;
    let Rvz;
    if(uy2uz2==0){
        Rvy=vx[1];
        Rvz=vx[2];
    }
    else{
        Rvy =(-vx[0]*ux[1]+(vx[1]*(ux[0]*ux[1]**2.0+Uabs*ux[2]**2.0)+vx[2]*ux[1]*ux[2]*(ux[0]-Uabs))/uy2uz2)/Uabs;
        Rvz =(-vx[0]*ux[2]+(vx[2]*(ux[0]*ux[2]**2.0+Uabs*ux[1]**2.0)+vx[1]*ux[2]*ux[1]*(ux[0]-Uabs))/uy2uz2)/Uabs;
    }
    
    let cosphi=Rvy/Math.sqrt(Rvy**2.0+Rvz**2.0);
    let sinphi=Rvz/Math.sqrt(Rvy**2.0+Rvz**2.0);

    Rot[0][0]= ux[0]/Uabs;
    Rot[0][1]= ux[1]/Uabs;
    Rot[0][2]= ux[2]/Uabs;
    Rot[1][0]= -(ux[1]*cosphi+ux[2]*sinphi)/Uabs;
    Rot[1][1]=(ux[2]*(ux[2]*cosphi-ux[1]*sinphi)+ux[0]*ux[1]*(ux[1]*cosphi+ux[2]*sinphi)/Uabs)/uy2uz2;
    Rot[1][2]=(ux[1]*(ux[1]*sinphi-ux[2]*cosphi)+ux[0]*ux[2]*(ux[2]*sinphi+ux[1]*cosphi)/Uabs)/uy2uz2;
    Rot[2][0]=(ux[1]*sinphi-ux[2]*cosphi)/Uabs;
    Rot[2][1]=(-ux[2]*(ux[1]*cosphi+ux[2]*sinphi)+ux[0]*ux[1]*(ux[2]*cosphi-ux[1]*sinphi)/Uabs)/uy2uz2;
    Rot[2][2]=(ux[1]*(ux[1]*cosphi+ux[2]*sinphi)+ux[0]*ux[2]*(ux[2]*cosphi-ux[1]*sinphi)/Uabs)/uy2uz2;

    for (let i=0;i<3;i++){
        a_unit[i]= Rot[i][0];
        b_unit[i]= Rot[i][0]*Math.cos(gamma)+Rot[i][1]*Math.sin(gamma);
        c_unit[i]= Rot[i][0]*Math.cos(beta)+Rot[i][1]*DD+Rot[i][2]*PP;
    } 

    // output parameters: 3 reciprocal lattice vectors, a*, b*, and c*
    for (let i=0;i<3;i++){
        a_star[i]= 2.0*Math.PI/a/PP/Math.sin(gamma)*(b_unit[(i+1)%3]*c_unit[(i+2)%3]-b_unit[(i+2)%3]*c_unit[(i+1)%3]);
        b_star[i]= 2.0*Math.PI/b/PP/Math.sin(gamma)*(c_unit[(i+1)%3]*a_unit[(i+2)%3]-c_unit[(i+2)%3]*a_unit[(i+1)%3]);
        c_star[i]= 2.0*Math.PI/c/PP/Math.sin(gamma)*(a_unit[(i+1)%3]*b_unit[(i+2)%3]-a_unit[(i+2)%3]*b_unit[(i+1)%3]);
    }

}

function draw_DetMap(){

    var canvas = document.getElementById('CanvasDetMap');
    var context = canvas.getContext('2d');

    //refresh
    context.clearRect(0, 0, canvas.width, canvas.height);
    context.strokeStyle = "rgb(0, 0, 0)";
    context.lineWidth=1;

    //set background color
    context.fillStyle = "rgb(0, 0, 100)";
    context.fillRect(0, 0, canvas.width, canvas.height);

    Hmax = Number(document.getElementById('Hmax').value);
    Kmax = Number(document.getElementById('Kmax').value);
    Lmax = Number(document.getElementById('Lmax').value);

    // line
//    context.strokeStyle = "rgb(255, 0, 0)";
//    context.beginPath();
//    context.moveTo(X0, Y0);
//    context.lineTo(X0+200, Y0);
//    context.stroke();

    // color setting for circles indicating reflections
    context.strokeStyle = "rgb(250, 250, 0)";
    context.fillStyle = "rgb(250, 250, 0)";
    context.lineWidth=2;

    var Ghkl=new Array(3);
    for (var H=-Hmax;H<=Hmax;H+=1){
        for (var K=-Kmax;K<=Kmax;K+=1){
            for (var L=-Lmax;L<=Lmax;L+=1){

                if((H==0)&&(K==0)&&(L==0)){
                }
                else{
                    for(let i=0;i<3;i++){
                        Ghkl[i]=H*a_star[i]+K*b_star[i]+L*c_star[i];
                    }
    
                    if(Ghkl[0]==0.0){
                    }
                    else{
                        let G_sq = Ghkl[0]**2.0+Ghkl[1]**2.0+Ghkl[2]**2.0;
                        let Ki = Math.abs(0.5*G_sq/Ghkl[0]); // Ki >0
                        lambda = 2.0*Math.PI/Ki;
                        //lambda = Math.abs(4.0*Math.PI*Ghkl[0]/G_sq);

                        let phiv = Math.atan2(Ghkl[2], Math.sqrt((Ghkl[0]+Ki)**2.0+Ghkl[1]**2.0));
                        let phih = Math.atan2(Ghkl[1],Ghkl[0]+Ki);

                        let PosX=scaleX*phih/Math.PI*180.0/maxphih+X0;
                        let PosY=scaleY*(HD+LD/2-L20*Math.tan(phiv))/LD

                        context.beginPath();
                        context.arc(PosX,PosY, radius, 0, 2 * Math.PI);
                        context.stroke();
                    }  
                }
            }
        }
    }


//text for debug
//    context.font = "italic 13px sans-serif";
//    context.fillText(cosphi, X0, Y0);


}




function rot_Lattice(rot_ax_dir){
    let deg = 0.0;
    let xyz         // xyz=(0,1,2) for (x, y, z)-axis respectively.
    switch(rot_ax_dir){
        case 'rot_x_plus':
            deg = Number(document.getElementById('rot_x_deg').value)/180.0*Math.PI;
            xyz =0.0;
            break;
        case 'rot_x_minus':
            deg = (-1.0)*Number(document.getElementById('rot_x_deg').value)/180.0*Math.PI;
            xyz =0.0;
            break;
        case 'rot_y_plus':
            deg = Number(document.getElementById('rot_y_deg').value)/180.0*Math.PI;
            xyz =1.0;
            break;
        case 'rot_y_minus':
            deg = (-1.0)*Number(document.getElementById('rot_y_deg').value)/180.0*Math.PI;
            xyz =1.0;
            break;
        case 'rot_z_plus':
            deg = Number(document.getElementById('rot_z_deg').value)/180.0*Math.PI;
            xyz =2.0;
            break;
        case 'rot_z_minus':
            deg = (-1.0)*Number(document.getElementById('rot_z_deg').value)/180.0*Math.PI;
            xyz =2.0;
            break;
        default:
    }
    xyz_rotation(xyz,deg);
}

function xyz_rotation(xyz,deg){
    let r00;
    let r01;

    r00=a_star[(xyz+1)%3]*Math.cos(deg)-a_star[(xyz+2)%3]*Math.sin(deg);
    r01=a_star[(xyz+1)%3]*Math.sin(deg)+a_star[(xyz+2)%3]*Math.cos(deg);
    a_star[(xyz+1)%3]=r00;
    a_star[(xyz+2)%3]=r01;
    r00=b_star[(xyz+1)%3]*Math.cos(deg)-b_star[(xyz+2)%3]*Math.sin(deg);
    r01=b_star[(xyz+1)%3]*Math.sin(deg)+b_star[(xyz+2)%3]*Math.cos(deg);
    b_star[(xyz+1)%3]=r00;
    b_star[(xyz+2)%3]=r01;
    r00=c_star[(xyz+1)%3]*Math.cos(deg)-c_star[(xyz+2)%3]*Math.sin(deg);
    r01=c_star[(xyz+1)%3]*Math.sin(deg)+c_star[(xyz+2)%3]*Math.cos(deg);
    c_star[(xyz+1)%3]=r00;
    c_star[(xyz+2)%3]=r01;
}

function draw_OriViewer(){
    // サイズを指定
    const width = 800;
    const height = 400;
  
    // レンダラーを作成
    const renderer = new THREE.WebGLRenderer({
      canvas: document.querySelector('#OrientationViewer'),
      antialias: true
    });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(width, height);
    renderer.setClearColor(0xf8f8f8);
  
    // シーンを作成
    const scene = new THREE.Scene();
  
    // カメラを作成
    const camera = new THREE.PerspectiveCamera(30, width / height);
    camera.position.set(-800, 800, 800);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  
    //note: 
    //HRC coordinates (x(||ki),y,z(vertical)) 
    //THREE.js coordinates (x3,y3,z3) 
    //transformation : x3=x, y3=z, z3=-y  
    // detector bank 1
    const geometry1 = new THREE.BoxGeometry(50.0, LD*DetBankScale, DetBankWidth*DetBankScale);
    const material1 = new THREE.MeshStandardMaterial({ color: 0xC0C0C0 });  // color of detector bank
    const mesh1 = new THREE.Mesh(geometry1, material1);
    scene.add(mesh1);
    mesh1.rotation.y += DetBankAngles[0];       //rotation about the y axis.
    mesh1.position.x += L20*DetBankScale*Math.cos(DetBankAngles[0]);    // move along the x axis.
    mesh1.position.z -= L20*DetBankScale*Math.sin(DetBankAngles[0]);    // move along the z axis.
  
    // detector bank 2
    const geometry2 = new THREE.BoxGeometry(50.0, LD*DetBankScale, DetBankWidth*DetBankScale);
    const mesh2 = new THREE.Mesh(geometry2, material1);
    scene.add(mesh2);
    mesh2.rotation.y += DetBankAngles[1];
    mesh2.position.x += L20*DetBankScale*Math.cos(DetBankAngles[1]);
    mesh2.position.z -= L20*DetBankScale*Math.sin(DetBankAngles[1]);
  
    // detector bank 3
    const geometry3 = new THREE.BoxGeometry(50.0, LD*DetBankScale, DetBankWidth*DetBankScale);
    const mesh3 = new THREE.Mesh(geometry3, material1);
    scene.add(mesh3);
    mesh3.rotation.y += DetBankAngles[2];
    mesh3.position.x += L20*DetBankScale*Math.cos(DetBankAngles[2]);
    mesh3.position.z -= L20*DetBankScale*Math.sin(DetBankAngles[2]);
  
    // detector bank 4
    const geometry4 = new THREE.BoxGeometry(50.0, LD*DetBankScale, DetBankWidth*DetBankScale);
    const mesh4 = new THREE.Mesh(geometry4, material1);
    scene.add(mesh4);
    mesh4.rotation.y += DetBankAngles[3];
    mesh4.position.x += L20*DetBankScale*Math.cos(DetBankAngles[3]);
    mesh4.position.z -= L20*DetBankScale*Math.sin(DetBankAngles[3]);
  
    // detector bank 4
    const geometry5 = new THREE.BoxGeometry(2000,50,50);
    const mesh5 = new THREE.Mesh(geometry5, material1);
    scene.add(mesh5);
    mesh5.position.x -= 1300;
    
    //draw a*, b*, c*
    //a*
    var dir = new THREE.Vector3( a_star[0],a_star[2], -a_star[1] );
    var origin = new THREE.Vector3( 0, 0, 0 );
    var arrow_len = dir.length()*arrow_scale;
    var hex = 0xff0000;
    var arrowHelper = new THREE.ArrowHelper( dir.normalize(), origin, arrow_len, hex ,arrow_HeadLen,arrow_HeadWidth);
    scene.add(arrowHelper);
  
    //b*
    dir = new THREE.Vector3( b_star[0],b_star[2], -b_star[1] );
    arrow_len = dir.length()*arrow_scale;
    hex = 0x00ff00;
    arrowHelper = new THREE.ArrowHelper( dir.normalize(), origin, arrow_len, hex ,arrow_HeadLen,arrow_HeadWidth);
    scene.add(arrowHelper);
 
    //c*
    dir = new THREE.Vector3( c_star[0],c_star[2], -c_star[1] );
    arrow_len = dir.length()*arrow_scale;
    hex = 0x0000ff;
    arrowHelper = new THREE.ArrowHelper( dir.normalize(), origin,arrow_len, hex ,arrow_HeadLen,arrow_HeadWidth);
    scene.add(arrowHelper);
  
  
  
    // 平行光源
    const directionalLight = new THREE.DirectionalLight(0xffffff);
    directionalLight.position.set(-150, 240, 500);
    scene.add(directionalLight);
  
    const light = new THREE.AmbientLight(0xffffff, 1.0);
    scene.add(light);  
    // ポイント光源
  //  const pointLight = new THREE.PointLight(0xffffff, 2, 1000);
  //  scene.add(pointLight);
  //  const pointLightHelper = new THREE.PointLightHelper(pointLight, 3);
  //  scene.add(pointLightHelper);
  
    renderer.render(scene, camera);
  //  tick();
  
  }