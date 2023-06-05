public class FiniteElementMethod {

    class GlobalData {
        double Tot; // temperatura otoczenia
        double alfa; // warunek brzegowy konwekcji
        double q; // gęstość strumienia ciepła
        double k; // współczynnik przewodzenia ciepła
        double L; // długość pręta
        double S; // pole przekroju pręta
        int nN; // Liczba węzłów siatki MES
        int nE; // Liczba elementów siatki MES
    }

    class Node {
        double t; // temperatura na węźle
        double x; // obliczyć L zamiast x- dla elementu
        int BC; // warunki brzegowe 1 - konwekcja, 2 - stały strumień ciepła
    }

    class Element {
        double k; // współczynnik przewodzenia ciepła
        int[] id = new int[2]; // E1-N1, N2 - id[1.2]
        double Le; // długość elementu
        double[][] H; // macierz lokalna
        double[] P; // wektor lokalny

        public void initializeLocalMatrices() {
            H = new double[2][2];
            P = new double[2];
        }

        public void setLocalMatrixHValues(double value) {
            H[0][0] = value;
            H[0][1] = -value;
            H[1][0] = -value;
            H[1][1] = value;
        }

        public void setLocalVectorPValues(double value) {
            P[0] = value;
            P[1] = value;
        }
    }

    class SiatkaMES {
        Element[] E;
        Node[] N;

        SiatkaMES(int nE, int nN) {
            E = new Element[nE];
            N = new Node[nN];
            for (int i = 0; i < nE; i++) {
                E[i] = new Element();
            }
            for (int i = 0; i < nN; i++) {
                N[i] = new Node();
            }
        }
    }

    class SOE {
        double[][] HG; // nN x nN
        double[] PG; // nN x 1
        double[] T; // nN x 1
        int nN;

        SOE(int nN) {
            this.nN = nN;
        }

        public void initializeMatricesAndVectors(int nN) {
            HG = new double[nN][nN];
            PG = new double[nN];
        }

        public void addValuesToGlobalMatricesAndVectors(Element elem) {
            HG[elem.id[0] - 1][elem.id[0] - 1] += elem.H[0][0];
            HG[elem.id[0] - 1][elem.id[1] - 1] += elem.H[0][1];
            HG[elem.id[1] - 1][elem.id[0] - 1] += elem.H[1][0];
            HG[elem.id[1] - 1][elem.id[1] - 1] += elem.H[1][1];

            PG[elem.id[0] - 1] += elem.P[0];
            PG[elem.id[1] - 1] += elem.P[1];
        }

        void gaussElimination() {
            int n = HG.length;
            for (int p = 0; p < n; p++) {
                int max = p;
                for (int i = p + 1; i < n; i++) {
                    if (Math.abs(HG[i][p]) > Math.abs(HG[max][p])) {
                        max = i;
                    }
                }
                double[] temp = HG[p];
                HG[p] = HG[max];
                HG[max] = temp;
                double t = PG[p];
                PG[p] = PG[max];
                PG[max] = t;
                for (int i = p + 1; i < n; i++) {
                    double alpha = HG[i][p] / HG[p][p];
                    PG[i] -= alpha * PG[p];
                    for (int j = p; j < n; j++) {
                        HG[i][j] -= alpha * HG[p][j];
                    }
                }
            }
            T = new double[n];
            for (int i = n - 1; i >= 0; i--) {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++) {
                    sum += HG[i][j] * T[j];
                }
                T[i] = (PG[i] - sum) / HG[i][i];
            }
        }

        void handleBoundaryCondition(GlobalData data, SiatkaMES mesh, SOE equationSystem) {
            // Sprawdzenie warunków brzegowych dla każdego węzła
            for (int i = 0; i < data.nN; i++) {
                Node node = mesh.N[i];
                // Jeśli warunki brzegowe to konwekcja
                if (node.BC == 1) {
                    equationSystem.HG[i][i] += data.alfa * data.S;
                    equationSystem.PG[i] += data.alfa * data.Tot * data.S;
                }
                // Jeśli warunki brzegowe to stały strumień ciepła
                else if (node.BC == 2) {
                    // Przyjmujemy, że strumień ciepła jest dodawany do wektora P
                    equationSystem.PG[i] += data.q * data.S;
                }
            }
        }

        public void printResults() {
            System.out.println("Temperatury w węzłach siatki MES:");
            for (int i = 0; i < nN; i++) {
                System.out.println("Węzeł " + (i+1) + ": " + T[i]);
            }
        }
    }

    public void deriveTemperatures(GlobalData data, SiatkaMES grid, SOE systemOfEquations) {

        systemOfEquations.initializeMatricesAndVectors(data.nN);

        for (int i = 0; i < data.nE; i++) {
            Element elem = grid.E[i];
            Node node1 = grid.N[elem.id[0] - 1];
            Node node2 = grid.N[elem.id[1] - 1];

            elem.Le = node2.x - node1.x;
            elem.initializeLocalMatrices();

            double h1 = elem.k * data.S / elem.Le;
            elem.setLocalMatrixHValues(h1);

            double p1 = data.q * data.S * elem.Le / 2.0;
            elem.setLocalVectorPValues(p1);

            systemOfEquations.addValuesToGlobalMatricesAndVectors(elem);
        }

        systemOfEquations.handleBoundaryCondition(data, grid, systemOfEquations);

        systemOfEquations.gaussElimination();

        systemOfEquations.printResults();
    }


    public void run() {
        GlobalData inputData = new GlobalData();
        inputData.Tot = 400.0;
        inputData.alfa = 10.0;
        inputData.q = 150.0;
        inputData.k = 50.0;
        inputData.L = 5.0;
        inputData.S = 2.0;
        inputData.nN = 3;
        inputData.nE = 2;

        SiatkaMES network = new SiatkaMES(inputData.nE, inputData.nN);

        double dx = inputData.L / (inputData.nN - 1);
        for (int i = 0; i < inputData.nN; i++) {
            network.N[i].x = dx * i;
            network.N[i].t = 0.0;
            network.N[i].BC = 1;
        }

        for (int i = 0; i < inputData.nE; i++) {
            network.E[i].k = inputData.k;
            network.E[i].id[0] = i + 1;
            network.E[i].id[1] = i + 2;
        }

        SOE systemOfEquations = new SOE(inputData.nN);
        deriveTemperatures(inputData, network, systemOfEquations);
    }


}