import { Test, test } from 'tape'
import { SMILESInChIApi, InlineResponse200} from './typescript-node-client/api'

const basePath : string = "http://192.168.99.100:8080/v1"

function verifyInchi(inchi : String, t : Test) : void {
  t.equal("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", inchi)
}

function testSmilesToInchi(t : Test) : void {
    const api = new SMILESInChIApi(basePath)
    api.smilesToInchiGet("CCO")
      .then(value => verifyInchi(value.body.inchi, t))
      .error((error: any) => t.fail("Api call failed"))
}

test.test('smiles to inchi', testSmilesToInchi);
