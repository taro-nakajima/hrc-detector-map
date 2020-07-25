window.addEventListener('load', init);

function init() {
  // サイズを指定
  const width = 960;
  const height = 540;

  // レンダラーを作成
  const renderer = new THREE.WebGLRenderer({
    canvas: document.querySelector('#myCanvas'),
    antialias: true
  });
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setSize(width, height);
  renderer.setClearColor(0xf8f8f8);

  // シーンを作成
  const scene = new THREE.Scene();

  // カメラを作成
  const camera = new THREE.PerspectiveCamera(30, width / height);
  camera.position.set(0, 1200, +1100);
  camera.lookAt(new THREE.Vector3(0, 0, 0));

  // box
  const geometry = new THREE.BoxGeometry(300, 100, 64);
  // マテリアルを作成
  const material = new THREE.MeshStandardMaterial({ color: 0xC0C0C0 });
  // メッシュを作成
  const mesh = new THREE.Mesh(geometry, material);
  // 3D空間にメッシュを追加
  scene.add(mesh);


 mesh.rotation.y += 0.6;
//  mesh.rotation.x += 0.6;
  mesh.position.x += 100;
  //  mesh.rotation.z += 0.6;

  // 平行光源
  const directionalLight = new THREE.DirectionalLight(0xffffff);
  directionalLight.position.set(0, 100, 100);
  scene.add(directionalLight);

  const light = new THREE.AmbientLight(0xffffff, 1.0);
  scene.add(light);  
  // ポイント光源
//  const pointLight = new THREE.PointLight(0xffffff, 2, 1000);
//  scene.add(pointLight);
//  const pointLightHelper = new THREE.PointLightHelper(pointLight, 30);
//  scene.add(pointLightHelper);

  renderer.render(scene, camera);
//  tick();

  // 毎フレーム時に実行されるループイベントです
//  function tick() {
//    // メッシュを回転させる
//    mesh.rotation.x += 0.01;
//    mesh.rotation.y += 0.01;

    // ライトを周回させる
//    pointLight.position.set(
//      500 * Math.sin(Date.now() / 500),
//      500 * Math.sin(Date.now() / 1000),
//      500 * Math.cos(Date.now() / 500)
//    );

    // レンダリング
//    renderer.render(scene, camera);

//    requestAnimationFrame(tick);
//  }
}