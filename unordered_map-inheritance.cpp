

#include <iostream>
#include <valarray>
#include <cassert>
#include <algorithm>

#define XORRAND 3183507249

using std::unordered_map;
using std::function;
using std::cin;
using std::cout;
using std::string;

template<typename F, typename S>
class pair{
public:
    F first;
    S second;

    pair(F first_, S second_):first(first_), second(second_){}
};

template<typename E>
class dStack {
private:
    short int size;
    short int capacity;
    short int index;//for insert
    E **doubleStackArray;
public:
    dStack() : index(0), size(0), capacity(5) {
        //doubleStackArray = (E *) malloc(5 * sizeof(E));
        doubleStackArray = new E*[5];
    }

    dStack(const dStack<E>& target): size(target.size), capacity(target.capacity), index(target.index){
        doubleStackArray = new E*[capacity];
        for (int i = 0; i < target.size; ++i) {
            doubleStackArray[i]=target.doubleStackArray[i];
        }
    }

    ~dStack() {
        free(doubleStackArray);
        doubleStackArray = nullptr;
    }

    void alloc(int newCapacity) {
        E** tmp=new E*[newCapacity];
        for (int i=0; i<capacity; i++){
            tmp[i]=doubleStackArray[i];
        }
        delete []doubleStackArray;
        doubleStackArray=tmp;
        capacity = static_cast<short>(newCapacity);
    }

    void push(E* newObj) {
        //fool security
        if (newObj) {
            if (size >= capacity) alloc(capacity * 2);
            doubleStackArray[index++] = newObj;
            size++;
        }
    }

    E* &peekFirst() {
        return doubleStackArray[index - 2];
    }

    E* &peekLast() {
        return doubleStackArray[index - 1];
    }

    void popLast() {
        size = size - 1 > 0 ? size - 1 : 0;
        index = index - 1 > 0 ? index - 1 : 0;
    }

    bool empty() {
        //return index < 1;
        return size==0;
    }

    bool hasOne() {
        return size==1;
        //return (index == 1 && doubleStackArray);
    }
};

class abstract_data_t {
public:
    struct SpecialNode {
        const int m_key;
        int m_value;
        size_t m_hash;
        SpecialNode *next;
        SpecialNode *left;
        SpecialNode *right;

        explicit SpecialNode(const int k, int v, size_t h)
                : m_key(k), m_value(v), m_hash(h), next(nullptr), left(nullptr), right(nullptr) {
        }
    };
    virtual ~abstract_data_t() = 0;
    [[nodiscard]] virtual size_t hashEvaluate(int key) const = 0;
    virtual void put(int key, int value) = 0;
    virtual bool empty() = 0;
    [[nodiscard]] virtual int at(const int &key) const = 0;
    virtual size_t length() = 0;
    virtual size_t capacity() = 0;
    virtual void putAtRearrange(SpecialNode *target, SpecialNode **nha, int mem, SpecialNode **newList) = 0;
    virtual void order(SpecialNode *node, int newMem, SpecialNode **nha, SpecialNode **nlon) = 0;
    virtual void add(SpecialNode *&root, int key, int value, size_t hash) = 0;
    virtual void copyNode(SpecialNode *initOb, SpecialNode *where) = 0;
    virtual void haresize(int mem) = 0;
    virtual void clearRoot(SpecialNode *root) = 0;
    virtual void clearList(SpecialNode *nodeList) = 0;
    virtual void clear() = 0;
    virtual bool contains(int key) = 0;
    virtual size_t count(int key) = 0;
    virtual int &operator[](int key)=0;
    virtual bool operator==(abstract_data_t& b)=0;
    virtual bool operator!=(abstract_data_t& b)= 0;
};

inline abstract_data_t::~abstract_data_t() = default;


//template <typename K, typename V,  typename Iter = iterator<K,V>>
class UnorderedMap : public abstract_data_t{
private:
    int size;
    int m_capacity;
    const int load_factor = 1;

    class iterator {
    private:
        SpecialNode *currentGeneral;
        SpecialNode *currentLocal;
        dStack<SpecialNode> dobStack;
    public:
        explicit iterator(SpecialNode *&list) : currentGeneral(list) {
            findLastOnBranch(dobStack, currentGeneral);
            currentLocal = dobStack.peekLast();
        }

        iterator(iterator &iter)= default;

        iterator() : currentLocal(nullptr), currentGeneral(nullptr), dobStack() {};

        ~iterator() = default;

        iterator& operator++() {
            if (dobStack.hasOne())
                dobStack.popLast();

            if (dobStack.empty()) {
                currentGeneral = currentGeneral->next;//FIXME: initial map null
                if (currentGeneral) {
                    findLastOnBranch(dobStack, currentGeneral);
                    currentLocal = dobStack.peekLast();
                } else {
                    currentLocal = nullptr;
                    return *this;
                }
            } else {
                if (dobStack.peekFirst()->right == dobStack.peekLast()) {
                    currentLocal = dobStack.peekFirst();
                    dobStack.popLast();
                } else {
                    SpecialNode *tmp = dobStack.peekFirst();
                    dobStack.popLast();
                    if(tmp->right)
                        findLastOnBranch(dobStack, tmp->right);
                    currentLocal = dobStack.peekLast();
                }
            }
            return *this;
        }

        iterator operator++(int) {
            iterator copy(*this);
            ++(*this);
            return copy;
        }

        // is searching from the parent node including it
        // if x is not NULL I call "findLastOnBranch" with x as second parameter
        // searching has snake form
        template<typename E>
        void findLastOnBranch(dStack<E> &stack, SpecialNode *node) {
            if (!node->right && !node->left) {
                stack.push(node);
                return;
            }
            if (node->left) {
                stack.push(node);
                return findLastOnBranch(stack, node->left);
            }
            if (node->right) {
                stack.push(node);
                return findLastOnBranch(stack, node->right);
            }
        }

        pair<int, int> operator*() {
            return {currentLocal->m_key, currentLocal->m_value};
        }

        pair<int, int>* operator->(){
            return new pair(currentLocal->m_key, currentLocal->m_value);
        }

        bool operator!=(const iterator &other) const {
            return !(*this==other);
        }

        bool operator==(const iterator &other) const {
            return currentLocal == other.currentLocal;
        }
    };

    SpecialNode **hashArray;
    SpecialNode *listOfNodes;
public:

    iterator &begin() {
        iterator *tmp = new iterator(listOfNodes);
        return *tmp;
    };

    iterator &end() {
        return *new iterator();
    };

    iterator find(int key){
        for(auto i=this->begin(); i!=this->end(); i++){
            if(i->first==key)
                return i;
        }
        return this->end();
    }

    static bool findByKey(iterator begin, const iterator& end, int key){
        for(auto i=begin; i!=end; i++){
            if(i->first==key)
                return true;
        }
        return false;
    }

    UnorderedMap() : size(0), m_capacity(10), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
    };

    explicit UnorderedMap(int capacity) : size(0), m_capacity(capacity), hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
    };

    UnorderedMap(const int arr[], int length) : size(0), m_capacity(length * 2), hashArray(nullptr),
                                                listOfNodes(nullptr) {
        haresize(m_capacity);//*2
        for (int i = 0; i < length; ++i) {
            put(i, arr[i]);
        }
    };

    UnorderedMap(const UnorderedMap &initObject) : size(initObject.size), m_capacity(initObject.m_capacity),
                                                   hashArray(nullptr), listOfNodes(nullptr) {
        haresize(m_capacity);
        for (int i = 0; i < m_capacity; i++) hashArray[i] = nullptr;
        SpecialNode **initArray = initObject.hashArray;
        for (int i = m_capacity - 1; i >= 0; i--) {
            if (initArray[i]) {
                SpecialNode *tmp = new SpecialNode(initArray[i]->m_key, initArray[i]->m_value, initArray[i]->m_hash);
                hashArray[i] = tmp;
                copyNode(initArray[i], hashArray[i]);
            }
        }
        function<void(SpecialNode *, SpecialNode *&)> fillList = [&fillList, this](SpecialNode *fromList,
                                                                                   SpecialNode *&listON) -> void {
            if (fromList->next) {
                fillList(fromList->next, listON);
            }
            hashArray[(fromList->m_hash) % m_capacity]->next = listON;
            listON = hashArray[(fromList->m_hash) % m_capacity];
        };
        fillList(initObject.listOfNodes, listOfNodes);
    };

    UnorderedMap(iterator begin, const iterator& end): UnorderedMap(){
        for(auto i=begin; i!=end; i++){
            this->put(i->first,i->second);
        }
    }

    ~UnorderedMap() override{
        clear();
    };

    //for integers
    [[nodiscard]] size_t hashEvaluate(const int key) const override{
        size_t result = floor(((int) sqrt(key) ^ XORRAND) * (10 - M_PI / (key / 2 + key)));
        if (result > 10000) {
            std::string numberStr = std::to_string(result);
            int del = (numberStr[5] + numberStr[4]) * 10 + (numberStr[4] + numberStr[1]);
            result = floor((double) result * del);
        }
        return result;
    };

    //insert value in map by the key
    void put(const int key, const int value) override{
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
    };

    static bool is_equal(const UnorderedMap &object, UnorderedMap &target) {//meant const
        if (object.size == target.size && object.m_capacity == target.m_capacity) {
            SpecialNode *l1 = object.listOfNodes;
            while (l1) {
                if (!is_node_equal(l1, target)) return false;
                l1 = l1->next;
            }
            return true;
        }
        return false;
    };

    static int is_node_equal(SpecialNode *object, UnorderedMap& target) {//meant const
        if (target.find(object->m_key) == target.end()) return 0;
        int x = 1, y = 1;
        if (object->left) {
            x = is_node_equal(object->left, target);
        }
        if (object->right) {
            y = is_node_equal(object->right, target);
        }
        return x * y;
    }

    bool empty() override{
        return (!listOfNodes);
    };

    //swap containers' content
    void swap(UnorderedMap &other) {
        std::swap(m_capacity, other.m_capacity);
        std::swap(size, other.size);
        std::swap(hashArray, other.hashArray);
        std::swap(listOfNodes, other.listOfNodes);
    };

    //search by the key
    int at(const int &key) const override{
        int index = hashEvaluate(key) % m_capacity;
        function<SpecialNode *(SpecialNode *, const int)> find = [&find](SpecialNode *root,
                                                                         const int &key) -> SpecialNode * {
            if (root == nullptr) return nullptr;
            if (root->m_key == key) return root;
            if (key < root->m_key)
                return find(root->left, key);
            else if (key > root->m_key)
                return find(root->right, key);
            return nullptr;
        };
        SpecialNode *tmp = find(hashArray[index], key);
        if(!tmp){
            throw std::out_of_range("Out of range");
        }
        return tmp->m_value;
    };//what does it return

    //how many inserted elements
    size_t length() override{
        return size;
    }

    //how many elements
    size_t capacity() override{
        return m_capacity;
    }

    //rearrangement of nodes in struct after hash array resize(changing hash array, list of nodes)
    void putAtRearrange(SpecialNode *target, SpecialNode **nha, int mem, SpecialNode **newList) override{
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

    //adding node in the right place in mini-tree after hash array resize
    void add(SpecialNode *&root, int key, int value, size_t hash) override{
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
    void copyNode(SpecialNode *initOb, SpecialNode *where) override{
        if (initOb->left != nullptr) {
            SpecialNode *objectForInsert = new SpecialNode(initOb->left->m_key, initOb->left->m_value,
                                                           initOb->left->m_hash);
            where->left = objectForInsert;
            copyNode(initOb->left, where->left);
        }
        if (initOb->right != nullptr) {
            SpecialNode *objectForInsert = new SpecialNode(initOb->right->m_key, initOb->right->m_value,
                                                           initOb->right->m_hash);
            where->right = objectForInsert;
            copyNode(initOb->right, where->right);
        }
    }

    //hash array resize
    void haresize(int mem) override{// mem should be > m_capacity
        struct SpecialNode **newHashArray = new SpecialNode *[mem];
        for (int i = 0; i < mem; i++) newHashArray[i] = nullptr;
        if (!hashArray) {
            hashArray = newHashArray;
            m_capacity = mem;
            return;
        }
        SpecialNode *newListOfNodes = nullptr;
        for (int i = 0; i < m_capacity; i++) {
            if (hashArray[i]) {
                order(hashArray[i], mem, newHashArray, &newListOfNodes);
            }
        }
        listOfNodes = newListOfNodes;
        m_capacity = mem;
        delete[]hashArray;
        hashArray = newHashArray;
    };

    void order(SpecialNode *node, int newMem, SpecialNode **nha, SpecialNode **nlon) override{
        if (!node)
            return;
        if (node->left) {
            order(node->left, newMem, nha, nlon);
            node->left = nullptr;
        }
        if (node->right) {
            order(node->right, newMem, nha, nlon);
            node->right = nullptr;
        }
        putAtRearrange(node, nha, newMem, nlon);
    }

    //clearing mini-tree
    void clearRoot(SpecialNode *root) override{
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
    void clearList(SpecialNode *nodeList) override{
        if (nodeList == nullptr) return;
        if (nodeList->next == nullptr) {
            clearRoot(nodeList);
        } else {
            clearList(nodeList->next);
            nodeList->next = nullptr;
            clearRoot(nodeList);
        }
    }

    void clear() override{
        clearList(listOfNodes);
        listOfNodes = nullptr;
        for (int i = 0; i < m_capacity; ++i) {
            if (hashArray[i]) {
                hashArray[i] = nullptr;
            }
        }
        delete[]hashArray;
        hashArray = nullptr;
        m_capacity = 0;
        size = 0;
    }

    bool contains(int key) override{
        at(key);
        return true;
    };

    size_t count(int key) override{
        at(key);
        return 1;
    };

    bool operator==(abstract_data_t& b) override{//meant const
        return UnorderedMap::is_equal(*this, *dynamic_cast<UnorderedMap*>(&b));
    }

    bool operator!=(abstract_data_t& b) override {
        return !(*this==b);
    }


    UnorderedMap &operator=(const UnorderedMap &initObj) {
        UnorderedMap tmp(initObj);
        swap(tmp);
        tmp.clear();
        return *this;
    };

    int &operator[](int key) override{
        int index = hashEvaluate(key) % m_capacity;
        function<SpecialNode *(SpecialNode *, const int)> find = [&find](SpecialNode *root,
                                                                         const int &key) -> SpecialNode * {
            if (root == nullptr) return nullptr;
            if (root->m_key == key) return root;
            if (key < root->m_key)
                return find(root->left, key);
            else if (key > root->m_key)
                return find(root->right, key);
            return nullptr;
        };
        SpecialNode *tmp = find(hashArray[index], key);
        if (tmp == nullptr) {
            put(key, 0);
            tmp = find(hashArray[index], key);
            return tmp->m_value;
        }
        return tmp->m_value;
    }

    friend std::ostream &operator<<(std::ostream &stream, const UnorderedMap &a) {
        SpecialNode *tmp = a.listOfNodes;
        int *keys = static_cast<int *>(calloc(a.size, sizeof(int)));
        int i = 0;
        function<void(SpecialNode *)> findKeys = [&findKeys, &keys, &i](SpecialNode *root) -> void {
            if (root == nullptr) return;
            if (root->left != nullptr)
                findKeys(root->left);
            keys[i++] = root->m_key;
            if (root->right != nullptr)
                findKeys(root->right);
        };
        //finding keys
        for (; tmp; tmp = tmp->next) {
            findKeys(tmp);
        }
        //sorting keys
        for (i = 1; i < a.size; i++) {
            int targ = keys[i];
            int j = i - 1;
            while (j >= 0 && keys[j] > targ) {
                keys[j + 1] = keys[j];
                j = j - 1;
            }
            keys[j + 1] = targ;
        }
        //output
        for (i = 0; i < a.size; i++) {
            stream << a.at(keys[i]) << "\n";
        }
        free(keys);
        return stream;
    }

    friend std::istream &operator>>(std::istream &stream, UnorderedMap &a) {
        int tmp, i = 0;
        while (stream >> tmp) {
            a.put(i, tmp);
            i++;
        }
        return stream;
    }

};

//bool operator==(const abstract_data_t& a, const abstract_data_t& b){
//    return unordered_Map::is_equal(*dynamic_cast<const unordered_Map*>(&a), *dynamic_cast<const unordered_Map*>(&b));
// }
//bool operator!=(const abstract_data_t& a, const unordered_Map& b){
//    return !(a==b);
//}

void print(const abstract_data_t& a){
    std::cout<<*dynamic_cast<const UnorderedMap*>(&a);
}

int main() {
    const int a[4] = {2, 3, 5, 7};
    abstract_data_t *v = new UnorderedMap(a, 4);
    print(*v);
    assert(!v->empty());
    v->put(12,11);
    assert(11 == v->at(12));
    const int b[3] = {13, 17, 19};
    abstract_data_t *w = new UnorderedMap();
    assert(w->empty());
    *(UnorderedMap*)w = *(UnorderedMap*)v;
    assert(*w == *v);
    (*w)[0] = 0;
    assert(0 == (*w)[0]);
}
