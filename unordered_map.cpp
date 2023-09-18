#include <iostream>
#include <valarray>
#include <cassert>
#include <algorithm>
#include <string>

#define XORRAND 3183507249

using std::unordered_map;
using std::function;
using std::cin;
using std::cout;

template<typename K, typename V>//what's the difference between class and typename
class abstract_data_t {
private:
    int size;
    int m_capacity;
    const int load_factor = 1;

    struct SpecialNode {
        const K m_key;
        V m_value;
        size_t m_hash;
        SpecialNode *next;
        SpecialNode *left;
        SpecialNode *right;

        explicit SpecialNode(const K k, V v, size_t h)
                : m_key(k),m_value(v),m_hash(h),next(nullptr),left(nullptr),right(nullptr){
        }
    };

    SpecialNode **hashArray;
    SpecialNode *listOfNodes;
public:
    abstract_data_t() : size(0), m_capacity(10), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
    };

    abstract_data_t(int capacity) : size(0), m_capacity(capacity), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
    };

    abstract_data_t(const int arr[], int length) : size(0), m_capacity(length*2), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);//*2
        for (int i = 0; i < length; ++i) {
            put(i, arr[i]);
        }
    };

    abstract_data_t(const abstract_data_t &initObject) : size(initObject.size), m_capacity(initObject.m_capacity), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
        for (int i = 0; i < m_capacity; i++) hashArray[i]= nullptr;
        SpecialNode** initArray=initObject.hashArray;
        for (int i = m_capacity-1; i >=0; i--) {
            if (initArray[i]) {
                SpecialNode* tmp=new SpecialNode(initArray[i]->m_key, initArray[i]->m_value,initArray[i]->m_hash);
                hashArray[i]=tmp;
                copyNode(initArray[i], hashArray[i]);
            }
        }
        function<void(SpecialNode *,SpecialNode *&)> fillList = [&fillList, this](SpecialNode *fromList, SpecialNode *&listON) -> void{
            if(fromList->next){
                fillList(fromList->next, listON);
            }
            hashArray[(fromList->m_hash)%m_capacity]->next=listON;
            listON=hashArray[(fromList->m_hash)%m_capacity];
        };
        fillList(initObject.listOfNodes, listOfNodes);
    };

    ~abstract_data_t() {
        clear(listOfNodes);
        listOfNodes = nullptr;
        for (int i = 0; i < m_capacity; ++i) {
            if (hashArray[i]) {
                hashArray[i] = nullptr;
            }
        }
        delete []hashArray;
        hashArray = nullptr;
        m_capacity = 0;
        size = 0;
    };

    size_t hashEvaluate(K key) {
        size_t result = floor(((int) sqrt(key) ^ XORRAND) * (10 - M_PI / (key / 2 + key)));
        if (result > 10000) {
            std::string numberStr = std::to_string(result);
            int del = (numberStr[5] + numberStr[4]) * 10 + (numberStr[4] + numberStr[1]);
            result = floor(result * del);
        }
        return result;
    };

    //insert value in map by the key
    bool put(const K key, const V value) {
        if (size / m_capacity >= load_factor) haresize(m_capacity * 2);
        size_t hash = hashEvaluate(key);
        int index = hash % m_capacity;
        if (hashArray[index] == nullptr) {
            SpecialNode *newEl = new SpecialNode(key, value, hash);
            size++;
            hashArray[index] = newEl;
            newEl->next = listOfNodes;
            listOfNodes = newEl;
        } else {
            struct SpecialNode *tmp = hashArray[index];
            add(tmp, key, value, hash);
            size++;
        }
        return true;
    };

    static bool is_equal(const abstract_data_t &object, const abstract_data_t &target) {};// ???

    void assign(size_t amount, int repl) {};// ???

    bool empty() {
        if (listOfNodes) return false;
        return true;
    };

    //swaping containers' content
    void swap(abstract_data_t &other) {
        std::swap(m_capacity,other.m_capacity);
        std::swap(size,other.size);
        std::swap(hashArray,other.hashArray);
        std::swap(listOfNodes,other.listOfNodes);
    };

    //search by the key
    V at(const K &key) {
        int index = hashEvaluate(key) % m_capacity;
        function<SpecialNode *(SpecialNode *,const K)> find = [&find](SpecialNode *root, const K &key) ->  SpecialNode *{
            if (root == nullptr) return nullptr;
            if (root->m_key == key) return root;
            if (key < root->m_key)
                return find(root->left, key);
            else if (key > root->m_key)
                return find(root->right, key);
            return nullptr;
        };
        SpecialNode *tmp = find(hashArray[index], key);
        return tmp == nullptr ? 0 : tmp->m_value;
    };//what does it return

    //rearrangement of nodes in struct after hash array resize(changing hash array, list of nodes)
    void putAtRearrange(SpecialNode *target, SpecialNode **nha, int mem, SpecialNode **newList) {
        size_t hash = target->m_hash;
        int index = hash % mem;
        if (nha[index] == nullptr) {
            nha[index] = target;
            target->next = *newList;
            *newList = target;
        } else {
            struct SpecialNode *tmp = nha[index];
            add(tmp, target->m_key, target->m_value, hash);
        }
    }

    //rearrangement of nodes in struct after hash array resize(erasure of old pointers)
    void rearrange(SpecialNode *ha, int mem, SpecialNode **nha, SpecialNode **newList) {
        if (ha->left != nullptr) {
            rearrange(ha->left, mem, nha, newList);
            SpecialNode *rearOb = ha->left;
            ha->left = nullptr;
            putAtRearrange(rearOb, nha, mem, newList);
        }
        if (ha->right != nullptr) {
            rearrange(ha->right, mem, nha, newList);
            SpecialNode *rearOb = ha->right;
            ha->right = nullptr;
            putAtRearrange(rearOb, nha, mem, newList);
        }
    }

    //adding node in the right place in mini-tree after hash array resize
    void add(SpecialNode *&root, K key, V value, size_t hash) {
        if (root == nullptr) {
            root = new SpecialNode(key, value, hash);
        } else if (key < root->m_key) {
            add(root->left, key, value, hash);
        } else if (key > root->m_key) {
            add(root->right, key, value, hash);
        } else if (key == root->m_key) {
            root->m_value = value;
        }
    };

    //copying mini-tree
    void copyNode(SpecialNode* initOb, SpecialNode* where){
        if (initOb->left != nullptr) {
            SpecialNode *objectForInsert = new SpecialNode(initOb->left->m_key,initOb->left->m_value,initOb->left->m_hash);
            where->left=objectForInsert;
            copyNode(initOb->left, where->left);
        }
        if (initOb->right != nullptr) {
            SpecialNode *objectForInsert = new SpecialNode(initOb->right->m_key,initOb->right->m_value,initOb->right->m_hash);
            where->right=objectForInsert;
            copyNode(initOb->right, where->right);
        }
    }

    //hash array resize
    void haresize(int mem) {// mem should be > m_capacity
        struct SpecialNode **newHashArray = new SpecialNode *[mem];
        for (int i = 0; i < mem; i++) newHashArray[i] = nullptr;
        if (!hashArray) {
            hashArray = newHashArray;
            m_capacity = mem;
            return;
        }
        SpecialNode *newListOfNodes= nullptr;
        for (int i = 0; i < m_capacity; i++) {
            if (hashArray[i]) {
                putAtRearrange(hashArray[i], newHashArray, mem, &newListOfNodes);
                rearrange(hashArray[i], mem, newHashArray, &newListOfNodes);
            }
        }
        listOfNodes = newListOfNodes;
        m_capacity = mem;
        delete []hashArray;
        hashArray = newHashArray;
    };

    //clearing mini-tree
    void clearRoot(SpecialNode *root) {
        if (root == nullptr) return;
        if (root->left != nullptr) {
            clearRoot(root->left);
            root->left = nullptr;
        }
        if (root->right != nullptr) {
            clearRoot(root->right);
            root->right = nullptr;
        }
        delete root;
    }

    //clearing tree nodes
    void clear(SpecialNode *nodeList) {
        if (nodeList == nullptr) return;
        if (nodeList->next == nullptr) {
            clearRoot(nodeList);
        } else {
            clear(nodeList->next);
            nodeList->next = nullptr;
            clearRoot(nodeList);
        }
    }

    void output() {
        SpecialNode *tmp = listOfNodes;
        function<void(SpecialNode *)> print = [&print](SpecialNode *root) -> void {
            if (root == nullptr) return;
            if (root->left != nullptr)
                print(root->left);
            cout << "key " << root->m_key << " - " << "value " << +root->m_value << "\n";
            if (root->right != nullptr)
                print(root->right);
        };
        while (tmp) {
            print(tmp);
            tmp = tmp->next;
        }
    }
};
int main(){}
