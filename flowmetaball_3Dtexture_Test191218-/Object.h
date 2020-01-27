﻿#pragma once
#include <array>

//
// 図形データ
//
class Object
{
  // 頂点配列オブジェクト名
  GLuint vao;

  // 頂点バッファオブジェクト名
  GLuint vbo;

  // 頂点の数
  const GLsizei vertexcount;

  // インデックスの頂点バッファオブジェクト
  GLuint ibo;

  // インデックスの数
  const GLsizei indexcount;

public:

  // 頂点属性
  struct Vertex
  {
    // 位置
    GLfloat position[3];

    // 法線
    GLfloat normal[3];

    // デフォルトコンストラクタ
    //   デフォルトコンストラクタは初期値を伴わない変数宣言の時に呼び出されます．
    //   コンストラクタの本体では何もしません．
    Vertex()
    {}

    // 引数で初期化を行うコンストラクタ
    //   頂点の位置と法線ベクトルの初期値を引数で指定するときに呼び出されます.
    //   コンストラクタの本体では何もしません．
    Vertex(float px, float py, float pz, float nx, float ny, float nz)
      : position{ px, py, pz }, normal{ nx, ny, nz }
    {}
  };

  // コンストラクタ
  //   vertexcount: 頂点の数
  //   vertex: 頂点属性を格納した配列
  //   indexcount: 頂点のインデックスの要素数
  //   index: 頂点のインデックスを格納した配列
  Object(GLsizeiptr vertexcount, const Vertex *vertex,
    GLsizeiptr indexcount, const GLuint *index)
    : vertexcount(static_cast<GLsizei>(vertexcount))
    , indexcount(static_cast<GLsizei>(indexcount))
  {
    // 頂点配列オブジェクト
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // 頂点バッファオブジェクト
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
      vertexcount * sizeof(Vertex), vertex, GL_STATIC_DRAW);

    // 結合されている頂点バッファオブジェクトを in 変数から参照できるようにする
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
      static_cast<char *>(0) + sizeof vertex->position);
    glEnableVertexAttribArray(1);

    // インデックスの頂点バッファオブジェクト
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
      indexcount * sizeof(GLuint), index, GL_STATIC_DRAW);
  }

  // デストラクタ
  virtual ~Object()
  {
    // 頂点配列オブジェクトを削除する
    glDeleteBuffers(1, &vao);

    // 頂点バッファオブジェクトを削除する
    glDeleteBuffers(1, &vbo);

    // インデックスの頂点バッファオブジェクトを削除する
    glDeleteBuffers(1, &ibo);
  }

private:

  // コピーコンストラクタによるコピー禁止
  Object(const Object &o);

  // 代入によるコピー禁止
  Object &operator=(const Object &o);

public:

  // 描画
  void draw() const
  {
    // 描画する頂点配列オブジェクトを指定する
    glBindVertexArray(vao);

    // 三角形で描画する
    glDrawElements(GL_TRIANGLES, indexcount, GL_UNSIGNED_INT, 0);
  }
};
