//
// Created by Henrique on 21/06/2019.
//

#ifndef RAYTRACERCHALLENGE_OBJ_PARSER_H
#define RAYTRACERCHALLENGE_OBJ_PARSER_H

#include <cstring>
#include <vector>
#include <cstdlib>
#include "structs.h"

long int getFaceIndex(char* text, char** out){
    auto index = strtol(text, out, 0);
    while(out[0][0] != ' ' && out[0][0] != '\n'){
        out[0] = out[0] + 1;
    }
    return index;
}

//TODO: Pode retornar triangle mesh ao invés de um group para ficar mais otimizado
Group parseOBJFile(const char* file){
    FILE* fp = fopen (file, "r");

    char line [128];
    char* aux = nullptr;

    auto vertices = std::vector<Tuple>();

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
                auto indices = std::vector<long int>();

                //Coloca na lista todos os índices da linha
                auto index = getFaceIndex(&line[2], &aux);
                if(index > 0 && index <= vertices.size()){
                    indices.push_back(index);
                    while(aux[0] != '\n'){
                        index = getFaceIndex(aux, &aux);
                        if(index > 0 && index <= vertices.size()){
                            indices.push_back(index);
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
                for(unsigned int i = 1; i < indices.size() - 1; i++){
                    g->insert(new Triangle(vertices[indices[0] - 1], vertices[indices[i] - 1], vertices[indices[i + 1] - 1]));
                }

                printf("");
            }
        }
    }

    return group;
}

#endif //RAYTRACERCHALLENGE_OBJ_PARSER_H
