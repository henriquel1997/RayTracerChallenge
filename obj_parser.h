//
// Created by Henrique on 21/06/2019.
//

#ifndef RAYTRACERCHALLENGE_OBJ_PARSER_H
#define RAYTRACERCHALLENGE_OBJ_PARSER_H

#include <cstring>
#include <vector>
#include <cstdlib>
#include "structs.h"

struct VerticeAndNormalIndex{
    long int vertice;
    long int normal;
};

VerticeAndNormalIndex getVerticeAndNormalIndex(char* text, char** out){
    long int vertice = strtol(text, out, 0);

    long int normal = 0;
    if(out[0][0] == '/'){
        out[0] = out[0] + 1;
        while(out[0][0] != '/'){
            out[0] = out[0] + 1;
        }
        out[0] = out[0] + 1;
        normal = strtol(out[0], out, 0);
    }

    while(out[0][0] != ' ' && out[0][0] != '\n'){
        out[0] = out[0] + 1;
    }

    return VerticeAndNormalIndex{vertice, normal};
}

//TODO: Pode retornar triangle mesh ao invés de um group para ficar mais otimizado
Group parseOBJFile(const char* file){
    FILE* fp = fopen (file, "r");

    char line [128];
    char* aux = nullptr;

    auto vertices = std::vector<Tuple>();
    auto normais = std::vector<Tuple>();

    auto group = Group();

    unsigned long long posMaisRecente = 0;

    struct GroupAndPos{
        const char* name;
        unsigned long long pos;
    };

    auto listaGrupos = std::vector<GroupAndPos>();

    while(fgets(line, sizeof(line), fp) != nullptr){
        //Verifica se tem um espaço depois do primeiro caractere
        if(line[1] == ' '){
            //Verifica qual é o primeiro caractere
            if(line[0] == 'v'){
                //Vértice
                auto v1 = strtod(&line[2], &aux);
                auto v2 = strtod(aux, &aux);
                auto v3 = strtod(aux, nullptr);

                vertices.push_back(point(v1, v2, v3));
            }else if(line[0] == 'g'){
                //Grupo
                auto tamanho = strlen(line);
                auto nome = (char*) malloc(sizeof(char) * (tamanho - 3));
                strncpy(nome, line + 2, tamanho - 3);
                nome[tamanho - 3] = '\0';

                bool achou = false;
                for(auto groupAndPos: listaGrupos){
                    if(strcmp(groupAndPos.name, nome) == 0){
                        posMaisRecente = groupAndPos.pos;
                        achou = true;
                        break;
                    }
                }

                if(achou){
                    continue;
                }

                //Criando um novo grupo
                posMaisRecente = group.size() + 1;
                group.insert(new Group());
                listaGrupos.push_back(GroupAndPos{ nome, posMaisRecente });

            }else if(line[0] == 'f'){
                //Face
                auto indices = std::vector<VerticeAndNormalIndex>();

                //Coloca na lista todos os índices da linha
                auto verticeAndNormal = getVerticeAndNormalIndex(&line[2], &aux);

                if(verticeAndNormal.vertice > 0 && verticeAndNormal.vertice <= vertices.size()){
                    indices.push_back(verticeAndNormal);
                    while(aux[0] != '\n'){
                        verticeAndNormal = getVerticeAndNormalIndex(aux, &aux);
                        if(verticeAndNormal.vertice > 0 && verticeAndNormal.vertice <= vertices.size()){
                            indices.push_back(verticeAndNormal);
                        }else{
                            break;
                        }
                    }
                }else{
                    continue;
                }

                Group* g;
                if(posMaisRecente == 0){
                    g = &group;
                }else{
                    g = dynamic_cast<Group*>(group.get(posMaisRecente - 1));
                }

                //Fan Triangulation
                auto v1 = indices[0].vertice;
                auto n1 = indices[0].normal;
                for(unsigned int i = 1; i < indices.size() - 1; i++){
                    auto v2 = indices[i].vertice;
                    auto n2 = indices[i].normal;
                    auto v3 = indices[i + 1].vertice;
                    auto n3 = indices[i + 1].normal;

                    if(n1 == 0 && n2 == 0 && n3 == 0){
                        //Triângulo Comum
                        g->insert(new Triangle(vertices[v1 - 1], vertices[v2 - 1], vertices[v3 - 1]));
                    }else{
                        //Smooth Triangle
                        g->insert(new Triangle(vertices[v1 - 1], vertices[v2 - 1], vertices[v3 - 1], normais[n1 - 1], normais[n2 - 1], normais[n3 - 1]));
                    }
                }
            }
        }else if(line[2] == ' '){
            //Normais dos Vértices
            if(line[0] == 'v' && line[1] == 'n'){
                auto n1 = strtod(&line[3], &aux);
                auto n2 = strtod(aux, &aux);
                auto n3 = strtod(aux, nullptr);

                normais.push_back(vector(n1, n2, n3));
            }
        }
    }

    return group;
}

#endif //RAYTRACERCHALLENGE_OBJ_PARSER_H
