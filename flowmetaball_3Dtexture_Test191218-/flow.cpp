//
// アプリケーション本体
//
#include "GgApplication.h"

// 粒子群オブジェクト
#include "Blob.h"

// 球オブジェクト
#include "Object.h"
#include "Uniform.h"



// 標準ライブラリ
#include <memory>
#include <random>

//
// 光源
//

// 光源の材質と位置
//GgVector ambient;   //! 光源強度の環境光成分.
//GgVector diffuse;   //! 光源強度の拡散反射光成分.
//GgVector specular;  //! 光源強度の鏡面反射光成分.
//GgVector position;  //! 光源の位置.
constexpr GgSimpleShader::Light lightData =
{
  { 0.1f, 0.1f, 0.1f, 1.0f },
  { 1.0f, 1.0f, 1.0f, 1.0f },
  { 0.0f, 0.0f, 0.0f, 0.0f },
  { 3.0f, 4.0f, 5.0f, 0.0f }
};


//
// カメラ
//

// カメラの位置
constexpr GLfloat cameraPosition[] = { 0.0f, 0.0f, 5.0f };

// 目標点の位置
constexpr GLfloat cameraTarget[] = { 0.0f, 0.0f, 0.0f };

// カメラの上方向のベクトル
constexpr GLfloat cameraUp[] = { 0.0f, 1.0f, 0.0f };

// 画角
constexpr GLfloat cameraFovy(1.0f);

// 前方面と後方面の位置
constexpr GLfloat cameraNear(3.0f), cameraFar(7.0f);

//
// 粒子群
//

// 生成する粒子群の数
const int bCount(8);

// 生成する粒子群の中心位置の範囲
const GLfloat bRange(0.8f);

// 一つの粒子群の粒子数
const int pCount(2000);

// 一つの粒子群の中心からの距離の平均
const GLfloat pMean(0.0f);

// 一つの粒子群の中心からの距離の標準偏差
const GLfloat pDeviation(0.3f);

// アニメーションの繰り返し間隔
const double interval(5.0);


//
// 粒子群の生成
//
//   paticles 粒子群の格納先
//   count 粒子群の粒子数
//   cx, cy, cz 粒子群の中心位置
//   rn メルセンヌツイスタ法による乱数
//   mean 粒子の粒子群の中心からの距離の平均値
//   deviation 粒子の粒子群の中心からの距離の標準偏差
//
void generateParticles(Particles &particles, int count,
  GLfloat cx, GLfloat cy, GLfloat cz,
  std::mt19937 &rn, GLfloat mean, GLfloat deviation)
{
  // 一様実数分布
  //   [0, 2) の値の範囲で等確率に実数を生成する
  std::uniform_real_distribution<GLfloat> uniform(0.0f, 2.0f);

  // 正規分布
  //   平均 mean、標準偏差 deviation で分布させる
  std::normal_distribution<GLfloat> normal(mean, deviation);

  // 格納先のメモリをあらかじめ確保しておく
  particles.reserve(particles.size() + count);

  // 原点中心に直径方向に正規分布する粒子群を発生する
  for (int i = 0; i < count; ++i)
  {
    // 緯度方向
    const GLfloat cp(uniform(rn) - 1.0f);
    const GLfloat sp(sqrt(1.0f - cp * cp));

    // 経度方向
    const GLfloat t(3.1415927f * uniform(rn));
    const GLfloat ct(cos(t)), st(sin(t));

    // 粒子の粒子群の中心からの距離 (半径)
    const GLfloat r(normal(rn));

    // 粒子を追加する
    particles.emplace_back(r * sp * ct + cx, r * sp * st + cy, r * cp + cz);

    //particles.emplace_back(cx, cy, cz, r * sp * ct, r * sp * st, r * cp);
  }
}

//
// アプリケーションの実行
//
void GgApplication::run()
{
  // ウィンドウを作成する
  Window window("Flow");

  //
  // 粒子群オブジェクトの作成
  //

  // 乱数の種に使うハードウェア乱数
  //std::random_device seed;

  // メルセンヌツイスタ法による乱数の系列を設定する
  //std::mt19937 rn(seed());
  std::mt19937 rn(54321);

  // 粒子群データの初期値
  Particles initial;

  // 一様実数分布
  //   [-1.0, 1.0) の値の範囲で等確率に実数を生成する
  std::uniform_real_distribution<GLfloat> center(-bRange, bRange);

  // 発生する粒子群の数だけ繰り返す
  for (int i = 0; i < bCount; ++i)
  {
    // 点の玉中心位置
    const GLfloat cx(center(rn)), cy(center(rn)), cz(center(rn));

    // 中心からの距離に対して密度が正規分布に従う点の玉を生成する
    generateParticles(initial, pCount, cx, cy, cz, rn, pMean, pDeviation);
  }

  // 粒子群オブジェクトを作成する
  std::unique_ptr<const Blob> blob(new Blob(initial));


  //
  // 球（物体）の描画
  //

  // プログラムオブジェクトを作成する
  const GLuint sphereShader(ggLoadShader("sphere.vert", "sphere.frag"));

  // uniform 変数の場所を取得する
  const GLint sphereModelviewLoc(glGetUniformLocation(sphereShader, "modelview"));
  const GLint sphereProjectionLoc(glGetUniformLocation(sphereShader, "projection"));
  const GLint sphereNormalMatrixLoc(glGetUniformLocation(sphereShader, "normalMatrix"));
  const GLint LposLoc(glGetUniformLocation(sphereShader, "Lpos"));
  const GLint LambLoc(glGetUniformLocation(sphereShader, "Lamb"));
  const GLint LdiffLoc(glGetUniformLocation(sphereShader, "Ldiff"));
  const GLint LspecLoc(glGetUniformLocation(sphereShader, "Lspec"));
  const GLint spherePointLoc(glGetUniformLocation(sphereShader, "position"));
  const GLint materialLoc(glGetUniformBlockIndex(sphereShader, "Material"));

  // uniform block の場所を 0 番の結合ポイントに結びつける
  glUniformBlockBinding(sphereShader, materialLoc, 0);

  // 球の半径
  constexpr float sphereRadius(0.2f);

  // 球の分割数
  const int slices(32), stacks(16);

  // 頂点属性を作る
  std::vector<Object::Vertex> solidSphereVertex;

  for (int j = 0; j <= stacks; ++j)
  {
    const float t(static_cast<float>(j) / static_cast<float>(stacks));
    const float y(sphereRadius * cos(3.141593f * t)), r(sphereRadius * sin(3.141593f * t));

    for (int i = 0; i <= slices; ++i)
    {
      const float s(static_cast<float>(i) / static_cast<float>(slices));
      const float z(r * cos(6.283185f * s)), x(r * sin(6.283185f * s));

      // 頂点属性を追加する
      solidSphereVertex.emplace_back(x, y, z, x, y, z);
    }
  }

  // インデックスを作る
  std::vector<GLuint> solidSphereIndex;

  for (int j = 0; j < stacks; ++j)
  {
    const int k((slices + 1) * j);

    for (int i = 0; i < slices; ++i)
    {
      // 頂点のインデックス
      const GLuint k0(k + i);
      const GLuint k1(k0 + 1);
      const GLuint k2(k1 + slices);
      const GLuint k3(k2 + 1);

      // 左下の三角形
      solidSphereIndex.emplace_back(k0);
      solidSphereIndex.emplace_back(k2);
      solidSphereIndex.emplace_back(k3);

      // 右上の三角形
      solidSphereIndex.emplace_back(k0);
      solidSphereIndex.emplace_back(k3);
      solidSphereIndex.emplace_back(k1);
    }
  }
  // 図形データを作成する
  std::unique_ptr<const Object> sphere(new Object(
    solidSphereVertex.size(), solidSphereVertex.data(),
    solidSphereIndex.size(), solidSphereIndex.data()));

  //
  // 描画の設定
  //

  // 図形描画用の光源
  const GgSimpleShader::LightBuffer light(lightData);
  //GgSimpleShader simpleShader("simple.vert", "simple.frag");

  //
  // 材質データ
  //
  struct Material
  {
    // 環境光の反射係数
    alignas(16) std::array<GLfloat, 3> ambient;
    // 拡散反射係数
    alignas(16) std::array<GLfloat, 3> diffuse;
    // 鏡面反射係数
    alignas(16) std::array<GLfloat, 3> specular;
    // 輝き係数
    alignas(4) GLfloat shininess;
  };


  // 光源データ
  static constexpr int Lcount(2);
  static constexpr Vector Lpos[] = { 0.0f, 0.0f, 5.0f, 1.0f, 8.0f, 0.0f, 0.0f, 1.0f };
  static constexpr GLfloat Lamb[] = { 0.2f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f };
  static constexpr GLfloat Ldiff[] = { 1.0f, 0.5f, 0.5f, 0.9f, 0.9f, 0.9f };
  static constexpr GLfloat Lspec[] = { 1.0f, 0.5f, 0.5f, 0.9f, 0.9f, 0.9f };

  // 色データ
  static constexpr Material color[] =
  {
    //      Kamb               Kdiff              Kspec        Kshi
    { 0.6f, 0.6f, 0.2f,  0.6f, 0.6f, 0.2f,  0.3f, 0.3f, 0.3f,  30.0f },
  { 0.1f, 0.1f, 0.5f,  0.1f, 0.1f, 0.5f,  0.4f, 0.4f, 0.4f,  60.0f }
  };
  const Uniform<Material> material(color, 2);

  // ビュー変換行列
  const GgMatrix mv(ggLookat(cameraPosition, cameraTarget, cameraUp));

  // 投影変換行列
  const GgMatrix mp(ggPerspective(cameraFovy, 1.0f, cameraNear, cameraFar));

  // 図形データ
  //const GgSimpleObj object("AC_1038.obj", simpleShader);
  //const GgSimpleObj object("bunny.obj", shader, true);
  //const GgSimpleObj object("box.obj", shader);

  // ビュー変換行列を求める
  const GgMatrix view(ggLookat(0.0f, 0.0f, 5.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));

  // 背面カリングを有効にする
  glEnable(GL_CULL_FACE);

  // デプスバッファを有効にする
  glEnable(GL_DEPTH_TEST);

  // 背景色を指定する
  glClearColor(0.1f, 0.2f, 0.3f, 0.0f);

  // 点のサイズはシェーダから変更する
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  /*
  ** レンダリング結果のブレンド
  **
  **   glBlendFunc(GL_ONE, GL_ZERO);                       // 上書き（デフォルト）
  **   glBlendFunc(GL_ZERO, GL_ONE);                       // 描かない
  **   glBlendFunc(GL_ONE, GL_ONE);                        // 加算
  **   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  // 通常のアルファブレンデング
  **   glBlendColor(0.01f, 0.01f, 0.01f, 0.0f);            // 加算する定数
  **   glBlendFunc(GL_CONSTANT_COLOR, GL_ONE);             // 定数を加算
  */

  // フレームバッファに加算する
  /**/
  // ポイントスプライトの設定
  // アルファブレンディングを設定する
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBlendFunc(GL_ONE, GL_ONE);
  glBlendEquation(GL_FUNC_ADD);

#if !DEPTH
  // デプスバッファは使わない
  glDisable(GL_DEPTH_TEST);
#endif

  // 時計をリセットする
  glfwSetTime(0.0);

  //
  // 描画
  //
  while (window)
  {
    // 定期的に粒子群オブジェクトをリセットする
    if (glfwGetTime() > interval)
    {
      blob->initialize(initial);
      glfwSetTime(0.0);
    }

    // ウィンドウを消去する
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // 透視投影変換行列を求める
    const GgMatrix projection(ggPerspective(1.0f, window.getAspect(), 1.0f, 10.0f));

    // モデル変換行列を求める
    const GgMatrix model(window.getTrackball());

    // モデルビュー変換行列を求める
    const GgMatrix modelview(view * model);

    // 粒子群オブジェクトを描画する
    //blob->draw(projection, modelview);

    //メタボールの描画
    blob->drawMetaball(projection, modelview, window, light);

    //
    // 球描画
    //

    //
    // 材質データ
    //
    struct Material
    {
      // 環境光の反射係数
      alignas(16) std::array<GLfloat, 3> ambient;

      // 拡散反射係数
      alignas(16) std::array<GLfloat, 3> diffuse;

      // 鏡面反射係数
      alignas(16) std::array<GLfloat, 3> specular;

      // 輝き係数
      alignas(4) GLfloat shininess;
    };

    // 色データ
    static constexpr Material color[] =
    {
      //      Kamb               Kdiff              Kspec        Kshi
      { 0.6f, 0.6f, 0.2f,  0.6f, 0.6f, 0.2f,  0.3f, 0.3f, 0.3f,  30.0f },
    { 0.1f, 0.1f, 0.5f,  0.1f, 0.1f, 0.5f,  0.4f, 0.4f, 0.4f,  60.0f }
    };
    const Uniform<Material> material(color, 2);

    //球の速度(倍速の値)
    const float Dspeed = 0.7f;

    // 初期位置
    float sphereX = 0.0f;
    float sphereY = 0.0f;
    float sphereZ = 5.0f - glfwGetTime() * Dspeed;

    //雲がちぎれる表現                    //球の中心座標 x        y        z
    const GgMatrix animation(modelview * ggTranslate(sphereX, sphereY, sphereZ));

    // 法線ベクトルの変換行列の格納先
    //GLfloat normalMatrix[9];

    // 法線ベクトルの変換行列を求める
    const GgMatrix normalMatrix(animation.normal());

    // オブジェクトのシェーダプログラムの使用開始
    glUseProgram(sphereShader);

    // uniform 変数に値を設定する
    glUniformMatrix4fv(sphereProjectionLoc, 1, GL_FALSE, projection.get());
    glUniformMatrix4fv(sphereModelviewLoc, 1, GL_FALSE, animation.get());
    //glUniformMatrix3fv(sphereNormalMatrixLoc, 1, GL_FALSE, normalMatrix);
    // 法線変換行列を設定する
    glUniformMatrix3fv(sphereNormalMatrixLoc, 1, GL_FALSE, normalMatrix.get());

    // アルファブレンディングを無効にして図形を描画する
    material.select(0);
    glDisable(GL_BLEND);
    sphere->draw();

    //オブジェクトの描画
    //simpleShader.use(mp, mv);
    //object.draw();

    // 粒子群オブジェクトを更新する
    blob->update();

    // カラーバッファを入れ替える
    window.swapBuffers();
  }
}
