import { Test, test } from 'tape'
import { CASInChIApi, CASInChIKeyApi, SMILESInChIApi, SMILESCASApi, InChIInChIKeyApi, SMILESInChIKeyApi } from './typescript-node-client/api'

const basePath : string = "http://192.168.99.100:8080/v1"

const ethanolInchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
const ethanolSmiles = "CCO"
const ethanolCas = "64-17-5"
const ethanolInchiKey = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

function verifyInchi(inchi : string, t : Test) : void {
  const expected = ethanolInchi
  t.equal(inchi, expected, expected)
}

function verifySmiles(smiles : string, t : Test) : void {
  const expected = ethanolSmiles
  t.equal(smiles, expected, expected)
}

function verifyInchiKey(inchi : string, t : Test) : void {
  const expected = ethanolInchiKey
  t.equal(inchi, expected, expected)
}

function verifyCas(smiles : string, t : Test) : void {
  const expected = ethanolCas
  t.equal(smiles, expected, expected)
}

function testSmilesToInchi(t : Test) : void {
    const api = new SMILESInChIApi(basePath)
    api.smilesToInchiGet(ethanolSmiles)
      .then(value => verifyInchi(value.body.inchi, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testInchiToSmiles(t : Test) : void {
    const api = new SMILESInChIApi(basePath)
    api.inchiToSmilesGet(ethanolInchi)
      .then(value => verifySmiles(value.body.smiles, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testSmilesToCas(t : Test) : void {
    const api = new SMILESCASApi(basePath)
    api.smilesToCasGet(ethanolSmiles)
      .then(value => verifyCas(value.body.cas, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testCasToSmiles(t : Test) : void {
    const api = new SMILESCASApi(basePath)
    api.casToSmilesGet(ethanolCas)
      .then(value => verifySmiles(value.body.smiles, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testInchiKeyToInchi(t : Test) : void {
    const api = new InChIInChIKeyApi(basePath)
    api.inchikeyToInchiGet(ethanolInchiKey)
      .then(value => verifyInchi(value.body.inchi, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testInchiToInchiKey(t : Test) : void {
    const api = new InChIInChIKeyApi(basePath)
    api.inchiToInchikeyGet(ethanolInchi)
      .then(value => verifyInchiKey(value.body.inchikey, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testSmilesToInchiKey(t : Test) : void {
    const api = new SMILESInChIKeyApi(basePath)
    api.smilesToInchikeyGet(ethanolSmiles)
      .then(value => verifyInchiKey(value.body.inchikey, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testInchiKeyToSmiles(t : Test) : void {
    const api = new SMILESInChIKeyApi(basePath)
    api.inchikeyToSmilesGet(ethanolInchiKey)
      .then(value => verifySmiles(value.body.smiles, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testCasToInchiKey(t : Test) : void {
    const api = new CASInChIKeyApi(basePath)
    api.casToInchikeyGet(ethanolCas)
      .then(value => verifyInchiKey(value.body.inchikey, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testInchiKeyToCas(t : Test) : void {
    const api = new CASInChIKeyApi(basePath)
    api.inchikeyToCasGet(ethanolInchiKey)
      .then(value => verifyCas(value.body.cas, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testCasToInchi(t : Test) : void {
    const api = new CASInChIApi(basePath)
    api.casToInchiGet(ethanolCas)
      .then(value => verifyInchi(value.body.inchi, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

function testInchiToCas(t : Test) : void {
    const api = new CASInChIApi(basePath)
    api.inchiToCasGet(ethanolInchi)
      .then(value => verifyCas(value.body.cas, t))
      .then(_ => t.end())
      .error((error: any) => t.fail("Api call failed"))
}

test.test('smiles to inchi', testSmilesToInchi);
test.test('inchi to smiles', testInchiToSmiles);
test.test('smiles to cas', testSmilesToCas);
test.test('cas to smiles', testCasToSmiles);
test.test('inchikey to inchi', testInchiKeyToInchi);
test.test('inchi to inchikey', testInchiToInchiKey);
test.test('smiles to inchikey', testSmilesToInchiKey);
test.test('inchikey to smiles', testInchiKeyToSmiles);
test.test('inchi to cas', testInchiToCas);
test.test('cas to inchi', testCasToInchi);
test.test('inchikey to cas', testInchiKeyToCas);
test.test('cas to inchikey', testCasToInchiKey);
